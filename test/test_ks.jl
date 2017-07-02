using Base.Test
using IMEXRKCB

# test for kuramoto-sivashinsky equation
Nₓ = 32
ν = (2π/39)^2

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

# initial condition
srand(0)
x₀ = 1e-2*randn(Nₓ)

# define scheme with quadrature
for scheme in [IMEXRK3R2R(IMEXRKCB3e, false, x₀),
               IMEXRK3R2R(IMEXRKCB3c, false, x₀),
               IMEXRK4R3R(IMEXRKCB4,  false, x₀)]

    # get integrator
    I = integrator(g, A, scheme, 1e-2)

    # warm up
    I(x₀, 10.0)

    # measure time and allocations
    t = @elapsed   I(x₀, 10.0)
    a = @allocated I(x₀, 10.0)

    # check
    @test t < 0.05
    @test a == 0
end