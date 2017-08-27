@testset "augmented state broadcast              " begin
    x = [1,   2,   3]
    q = [3.0, 4.0, 5.0]

    # define two augmented states
    z = IMEXRKCB.aug_state(x, q)
    y = IMEXRKCB.aug_state(copy(x)+1, copy(q)+1)

    # define some operation
    fun!(z, y) = (z .= 2.0.*z .+ 1.0.*y .- 1; z)

    # apply 
    fun!(z, y)

    # value
    @test IMEXRKCB._state(z) == [3, 6, 9]
    @test IMEXRKCB._quad(z)  == [9.0, 12.0, 15.0]

    # allocation
    @test (@allocated fun!(z, y)) == 0
end

@testset "quadrature                             " begin

    # define linear system
    g(t, x, ẋ) = (ẋ .= 0.5.*x; ẋ)
    A = Diagonal([0.5])

    # define example quadrature functions
    q(t, x, q̇) = (q̇[1] = 1; q̇[2] = x[1]; q̇[3] = t)

    # integration scheme
    for (impl, order, value) in [(IMEXRK3R2R(IMEXRKCB2,  false), 2, 19),
                                 (IMEXRK3R2R(IMEXRKCB3e, false), 3, 5.5),
                                 (IMEXRK3R2R(IMEXRKCB3c, false), 3, 5.8),
                                 (IMEXRK4R3R(IMEXRKCB4,  false), 4, 0.16)]

        # define scheme                                 
        scheme = IMEXRKScheme(impl, [0.0], [0.0, 0.0, 0.0])

        # exact values of integral
        exact = [5, exp(5) - exp(0), 5^2/2]

        # error should decrease at certain rate
        for Δt = logspace(0, -2.5, 5)

            # forward map
            f = integrator(g, A, q, scheme, Δt)

            # initial conditions
            x₀ = [1.0]
            q₀ = [0.0, 0.0, 0.0]
        
            # monitor the three quadrature equations
            mon = Monitor((xq->xq[2][1], 
                           xq->xq[2][2], 
                           xq->xq[2][3]), (x₀, q₀))

            # call
            f(x₀, q₀, 5, mon)

            # integrate 1 dt
            @test mon.samples[1] ≈ mon.times

            # integrate exp(t)dt
            @test norm(abs.(mon.samples[2] - (exp.(mon.times) - 1.0))) / Δt^(order-1) < value

            # integrate t*dt
            @test mon.samples[3] ≈ 0.5*(mon.times).^2

            # integrals
            @test norm(abs.(q₀ - exact)) / Δt^order < value
        end
    end
end