export IMEXMethod

struct IMEXMethod{X, TAB<:AbstractTableau, I}
    storage::Vector{X}
    tableau::TAB
    function IMEXMethod{I}(tableau::TAB, nstores::Int, x::X) where {X, TAB<:AbstractTableau}
        new{X, TAB, I}(X[similar(x) for i = 1:nstores], tableau)
    end
end

# allowed method tags
const __allowed__ = Dict{Symbol, Tuple{IMEXTableau, Int, Symbol}}(:CB2_3R2R  => (CB2,  3, :3R2R), 
                                                                  :CB3c_3R2R => (CB3c, 3, :3R2R), 
                                                                  :CB3e_3R2R => (CB3e, 3, :3R2R),
                                                                  :CB4_4R3R  => (CB4,  4, :4R3R))

# Main entry point
function IMEXMethod(tag::Symbol, x::X, q...)
    tag ∈ keys(__allowed__) || error("invalid method tag $tag")
    tableau, nstores, impl = __allowed__[tag]
    IMEXMethod{impl}(tableau, nstores, aug_state(x, q...))
end

# Three-register implementation of [2R] IMEXRK schemes from section 1.2.1 of CB 2015
function step!(method::IMEXMethod{X, TAB, :3R2R}, 
               sys::System, 
               t::Real, 
               Δt::Real, 
               x::X) where {X, TAB}

    # hoist temporaries out
    y   = method.storage[1]
    z   = method.storage[2]
    w   = method.storage[3]
    tab = method.tableau

    # loop over stages
    for k = 1:nstages(tab)
        if k == 1
            y .= x
        else
            y .= x .+ (tab[:aᴵ, k, k-1] - tab[:bᴵ, k-1])*Δt.*z 
                   .+ (tab[:aᴱ, k, k-1] - tab[:bᴱ, k-1])*Δt.*y
        end
        A_mul_B!(z, sys, y)                 # compute z = A*y then
        ImcA!(sys, tab[:aᴵ, k, k]*Δt, z, z) # get z = (I-cA)⁻¹*(A*y) in place
        w .= y .+ tab[:aᴵ, k, k]*Δt.*z      # w is the temp input, output is y
        sys(t + tab[:cᴱ, k]*Δt, w, y)
        x .= x .+ tab[:bᴵ, k]*Δt.*z .+ tab[:bᴱ, k]*Δt.*y
    end
    return nothing
end


# Four-register implementation of [3R] IMEXRK schemes from section 1.2.3 of CB 2015
function step!(method::IMEXMethod{X, TAB, :4R3R}, 
               sys::System, 
               t::Real, 
               Δt::Real, 
               x::X) where {X, TAB}

    # hoist temporaries out
    y   = method.store[1]
    zᴵ  = method.store[2]
    zᴱ  = method.store[3]
    w   = method.store[4]
    tab = method.tableau

    # loop over stages
    for k = 1:nstages(tab)
        if k == 1
             y .= x
            zᴵ .= x
            zᴱ .= x
        else
            zᴱ .= y .+ tab[:aᴱ, k, k-1]*Δt.*zᴱ
            if k < nstages(tab)
                y .= x .+ (tab[:aᴵ, k+1, k-1] .- tab[:bᴵ, k-1]).*Δt.*zᴵ 
                       .+ (tab[:aᴱ, k+1, k-1] .- tab[:bᴱ, k-1])./tab[:aᴱ, k+1, k-1].*(zᴱ .- y)
            end
            zᴱ .= zᴱ .+ tab[:aᴵ, k, k-1]*Δt.*zᴵ
        end
        A_mul_B!(w, sys, zᴱ)                 # compute w = A*zᴱ then
        ImcA!(sys, tab[:aᴵ, k, k]*Δt, w, zᴵ) # compute z = (I-cA)⁻¹*w in place
        w .= zᴱ .+ tab[:aᴵ, k, k]*Δt.*zᴵ     # w is the temp input, output is zᴱ
        sys(t + tab[:cᴱ, k]*Δt, w, zᴱ)
        x .= x .+ tab[:bᴵ, k]*Δt.*zᴵ .+ tab[:bᴱ, k]*Δt.*zᴱ        
    end
    return nothing
end