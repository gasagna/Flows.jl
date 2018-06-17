@testset "augmented state broadcast              " begin
    x = [1,   2,   3]
    q = [3.0, 4.0, 5.0]

    # define two augmented states
    z = Flows.aug_state(x, q)
    y = Flows.aug_state(copy(x)+1, copy(q)+1)

    # define some operation
    fun!(z, y) = (z .= 2.0.*z .+ 1.0.*y .- 1; z)

    # apply 
    fun!(z, y)

    # value
    @test Flows._state(z) == [3, 6, 9]
    @test Flows._quad(z)  == [9.0, 12.0, 15.0]

    # allocation
    @test (@allocated fun!(z, y)) == 0
end

@testset "quadrature                             " begin

    # define linear system
    g(t, x, ẋ) = (ẋ .= 0.5.*x; ẋ)
    A = Diagonal([0.5])

    # define example quadrature functions
    @inline quad(t, x, q̇) = (q̇[1] = 1; q̇[2] = x[1]; q̇[3] = t; return q̇)

    # state and quadrature
    x, q = Float64[0.0], Float64[0.0, 0.0, 0.0]

    # integration method
    for (method, order, value) in [(Scheme(:CB2_3R2R,  x, q), 2, 19),
                                   (Scheme(:CB3e_3R2R, x, q), 3, 5.5),
                                   (Scheme(:CB3c_3R2R, x, q), 3, 5.8),
                                   (Scheme(:CB4_4R3R,  x, q), 4, 0.16)]

        # exact values of integral
        exact = [5, exp(5) - exp(0), 5^2/2]

        # error should decrease at certain rate
        for Δt = logspace(0, -2.5, 5)

            # forward map
            f = integrator(g, A, quad, method, Δt)

            # initial conditions
            x₀ = Float64[1.0]
            q₀ = Float64[0.0, 0.0, 0.0]

            # call
            f(x₀, q₀, (0, 5))

            # integrals
            @test norm(abs.(q₀ - exact)) / Δt^order < value

            # test allocations
            fun(f, x₀, q₀, span) = @allocated f(x₀, q₀, span)
            @test fun(f, x₀, q₀, (0.0, 100.0)) <= 32
        end
    end
end