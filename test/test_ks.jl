using Base.Test
using IMEXRK

# test for kuramoto-sivashinsky

const Nₓ = 32
const ν = (2π/39)^2

# linear term
A = Diagonal(Float64[k*k*(1-ν*k*k) for k = 1:Nₓ])

# nonlinear term
function g(t::Real, x::AbstractVector, ẋ::AbstractVector)
    Nₓ = length(x)
    for k = 1:Nₓ
        s = zero(eltype(x))
        for m = max(-Nₓ, k-Nₓ):min(Nₓ, k+Nₓ)
            if !(k-m == 0 || m == 0)
                s += x[abs(m)]*x[abs(k-m)]*sign(m)*sign(k-m)
            end
        end
        ẋ[k] -= k*s
    end
    ẋ
end

# define system
f = DiagonalIMEXSystem(g, A)

# initial condition
srand(0)
x₀ = 1e-2*randn(Nₓ)

# define scheme
scheme = IMEXRK3R2R(IMEXRKCB3e, x₀)

# get flow
for Δt = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
    ϕ = forwmap!(f, 10, Δt, scheme)
    x = ϕ(copy(x₀))
    @show sum(abs2, x)
end