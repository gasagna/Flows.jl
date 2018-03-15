export IMEXMethod

# The third parameter is a symbol defining the implementation of the method,
# e.g. :3R2R for the 3 register implementation of 2R methods.
struct IMEXMethod{X, T<:Real, I}
    storage::Vector{X}
    tableau::IMEXTableau{T}
    function IMEXMethod{X, T, I}(tab::IMEXTableau{T}, nstores::Int, x::X) where {X, T, I}
        new(X[similar(x) for i = 1:nstores], tab)
    end
end

IMEXMethod(impl::Symbol, tab::IMEXTableau{T}, nstores::Int, x::X) where {X, T} =
    IMEXMethod{X, T, impl}(tab, nstores, x)

# allowed method tags. We convert tableaux to Float64. This should cover most cases.
D = Dict{Symbol, Tuple{IMEXTableau{Float64}, Int, Symbol}}
const __allowed__ =  D(:CB2_3R2R  => (convert(IMEXTableau{Float64}, CB2),  3, :_3R2R),
                       :CB3c_3R2R => (convert(IMEXTableau{Float64}, CB3c), 3, :_3R2R),
                       :CB3e_3R2R => (convert(IMEXTableau{Float64}, CB3e), 3, :_3R2R),
                       :CB4_4R3R  => (convert(IMEXTableau{Float64}, CB4),  4, :_4R3R))

# Main entry point
function IMEXMethod(tag::Symbol, x::X, q...) where {X}
    tag ∈ keys(__allowed__) || error("invalid method tag $tag")
    tab, nstores, impl = __allowed__[tag]
    IMEXMethod(impl, tab, nstores, aug_state(x, q...))
end

# Three-register implementation of [2R] IMEXRK schemes from section 1.2.1 of CB 2015
function step!(method::IMEXMethod{X, T, :_3R2R},
               sys::System,
               t::Real,
               Δt::Real,
               x::X) where {X, T}

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
            y .= x .+ (tab[:aᴵ, k, k-1] .- tab[:bᴵ, k-1]).*Δt.*z .+ 
                      (tab[:aᴱ, k, k-1] .- tab[:bᴱ, k-1]).*Δt.*y
        end
        A_mul_B!(z, sys, y)                 # compute z = A*y then
        ImcA!(sys, tab[:aᴵ, k, k]*Δt, z, z) # get z = (I-cA)⁻¹*(A*y) in place
        w .= y .+ tab[:aᴵ, k, k].*Δt.*z      # w is the temp input, output is y
        sys(t + tab[:cᴱ, k]*Δt, w, y)
        x .= x .+ tab[:bᴵ, k].*Δt.*z .+ tab[:bᴱ, k].*Δt.*y
    end
    return nothing
end


# Four-register implementation of [3R] IMEXRK schemes from section 1.2.3 of CB 2015
function step!(method::IMEXMethod{X, T, :_4R3R},
               sys::System,
               t::Real,
               Δt::Real,
               x::X) where {X, T}

    # hoist temporaries out
    y   = method.storage[1]
    zᴵ  = method.storage[2]
    zᴱ  = method.storage[3]
    w   = method.storage[4]
    tab = method.tableau

    # loop over stages
    for k = 1:nstages(tab)
        if k == 1
             y .= x
            zᴵ .= x
            zᴱ .= x
        else
            zᴱ .= y .+ tab[:aᴱ, k, k-1].*Δt.*zᴱ
            if k < nstages(tab)
                y .= x .+ (tab[:aᴵ, k+1, k-1] .- tab[:bᴵ, k-1]).*Δt.*zᴵ .+
                          (tab[:aᴱ, k+1, k-1] .- tab[:bᴱ, k-1])./tab[:aᴱ, k, k-1].*(zᴱ .- y)
            end
            zᴱ .= zᴱ .+ tab[:aᴵ, k, k-1]*Δt.*zᴵ
        end
        A_mul_B!(w, sys, zᴱ)                 # compute w = A*zᴱ then
        ImcA!(sys, tab[:aᴵ, k, k]*Δt, w, zᴵ) # compute z = (I-cA)⁻¹*w in place
        w .= zᴱ .+ tab[:aᴵ, k, k]*Δt.*zᴵ     # w is the temp input, output is zᴱ
        sys(t + tab[:cᴱ, k]*Δt, w, zᴱ)
        x .= x .+ tab[:bᴵ, k].*Δt.*zᴵ .+ tab[:bᴱ, k].*Δt.*zᴱ
    end
    return nothing
end