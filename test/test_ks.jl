# test for kuramoto-sivashinsky equation
const Nₓ = 32
const ν = (2π/39)^2

# linear term
const A = Diagonal(Float64[k*k*(1-ν*k*k) for k = 1:Nₓ])

# nonlinear term (the wrong way of doing it)
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

@testset "ks                                     " begin

    # initial condition
    srand(0)
    x₀ = 1e-2*randn(Nₓ)

    # define scheme with quadrature
    for (scheme, tmin) in [(Scheme(:CB2_3R2R,  x₀), 0.0195),
                           (Scheme(:CB3e_3R2R, x₀), 0.0262),
                           (Scheme(:CB3c_3R2R, x₀), 0.0261),
                           (Scheme(:CB4_4R3R,  x₀), 0.0410)]

        # get integrator
        I = integrator(g, A, scheme, 1e-2)

        # warm up
        I(x₀, (0.0, 10.0))

        # measure time and allocations
        t = minimum([@elapsed I(x₀, (0, 10)) for i = 1:100])
        a = @allocated I(x₀, (0, 10))

        # check
        @test t < tmin
        @test a == 0
    end
end