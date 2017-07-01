using Base.Test
using IMEXRKCB

@testset "quadrature                             " begin

    # define linear system
    g(t, x, ẋ) = (ẋ .= 0.5*x; ẋ)
    A = Diagonal([0.5])

    # define example quadrature function
    q(t, x, q̇) = (q̇[1] = 1; q̇[2] = x[1]; q̇[3] = t)

    # integration scheme
    for (scheme, order, value) in [(IMEXRK3R2R(IMEXRKCB3e, false, [0.0, 0.0], [0.0, 0.0, 0.0]), 3, 5.5),
                                   (IMEXRK3R2R(IMEXRKCB3c, false, [0.0, 0.0], [0.0, 0.0, 0.0]), 3, 5.8),
                                   (IMEXRK4R3R(IMEXRKCB4,  false, [0.0, 0.0], [0.0, 0.0, 0.0]), 4, 0.16)]

        # exact values of integral
        exact = [5, exp(5) - exp(0), 5^2/2]

        # error should decrease at certain rate
        for Δt = logspace(0, -2.5, 5)

            # forward map
            f = integrator(g, A, q, scheme, Δt)

            # initial conditions
            x₀ = [1.0]
            q₀ = [0.0, 0.0, 0.0]
        
            # monitor the quadrature equations
            ms = (Monitor((x, q)->q[1], [0.0], [0.0]),
                  Monitor((x, q)->q[2], [0.0], [0.0]),
                  Monitor((x, q)->q[3], [0.0], [0.0]))

            # call
            f(x₀, q₀, 5, ms...)

            # integrate 1 dt
            @test ms[1].samples ≈ ms[1].times

            # integrate exp(t)dt
            @test norm(abs.(ms[2].samples - (exp.(ms[2].times) - 1.0))) / Δt^(order-1) < value

            # integrate t*dt
            @test ms[3].samples ≈ 0.5*(ms[3].times).^2

            # integrals
            @test norm(abs.(q₀ - exact)) / Δt^order < value
        end
    end
end