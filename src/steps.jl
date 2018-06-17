export Scheme

# The third parameter is a symbol defining the implementation of the scheme,
# e.g. :3R2R for the 3 register implementation of 2R schemes.
struct Scheme{X, I, T<:Real, TAB<:AbstractTableau{T}}
    storage::Vector{X}
    tableau::TAB
    Scheme{I}(tab::TAB, nstores::Int, x::X) where {X, I, T, TAB<:AbstractTableau{T}} =
        new{X, I, T, TAB}(X[similar(x) for i = 1:nstores], tab)
end

# allowed scheme tags. We convert tableaux to Float64. This should cover most cases.
D = Dict{Symbol, Tuple{AbstractTableau, Int, Symbol}}
const __allowed__ =  D(:CB2_3R2R  => (convert(IMEXTableau{Float64}, CB2),  3, :_3R2R),
                       :CB3c_3R2R => (convert(IMEXTableau{Float64}, CB3c), 3, :_3R2R),
                       :CB3e_3R2R => (convert(IMEXTableau{Float64}, CB3e), 3, :_3R2R),
                       :CB4_4R3R  => (convert(IMEXTableau{Float64}, CB4),  4, :_4R3R))

# Main entry point
function Scheme(tag::Symbol, x::X, q...) where {X}
    tag ∈ keys(__allowed__) || error("invalid scheme tag $tag")
    tab, nstores, impl = __allowed__[tag]
    Scheme{impl}(tab, nstores, aug_state(x, q...))
end

# Three-register implementation of [2R] IMEXRK schemes from section 1.2.1 of CB 2015
function step!(scheme::Scheme{X, :_3R2R},
               sys::System,
               t::Real,
               Δt::Real,
               x::X) where {X}

    # hoist temporaries out
    y   = scheme.storage[1]
    z   = scheme.storage[2]
    w   = scheme.storage[3]
    tab = scheme.tableau

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
function step!(scheme::Scheme{X, :_4R3R},
               sys::System,
               t::Real,
               Δt::Real,
               x::X) where {X}

    # hoist temporaries out
    y   = scheme.storage[1]
    zᴵ  = scheme.storage[2]
    zᴱ  = scheme.storage[3]
    w   = scheme.storage[4]
    tab = scheme.tableau

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