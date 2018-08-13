@testset "coupled state broadcast                " begin
    x = [1,   2,   3]
    q = [3.0, 4.0, 5.0]

    # define two coupled states
    z = Flows.couple(x, q)
    y = Flows.couple(copy(x)+1, copy(q)+1)

    # define some operation
    fun!(z, y) = (z .= 2.0.*z .+ 1.0.*y .- 1; z)

    # apply 
    fun!(z, y)

    # value
    @test first(z) == [3, 6, 9]
    @test  last(z) == [9.0, 12.0, 15.0]

    # allocation
    @test (@allocated fun!(z, y)) == 0
end

@testset "quadrature                             " begin

    # define linear system
    g(t, x, ẋ) = (ẋ .= 0.5.*x; ẋ)
    A = Diagonal([0.5])

    # also define full explicit term for RK4
    gfull(t, x, ẋ) = (ẋ .= x; ẋ)

    # define example quadrature functions
    @inline quad(t, x, q̇) = (q̇[1] = 1; q̇[2] = x[1]; q̇[3] = t; return q̇)

    # state and quadrature
    x, q = Float64[0.0], Float64[0.0, 0.0, 0.0]

    # monitor first quadrature equation
    mon = Monitor(couple(x, q), copy)

    # integration method
    for (method, order, value, _g, _A) in [(RK4(     x, q, :NORMAL), 4, 6.2, gfull, nothing),
                                           (CB3R2R2( x, q, :NORMAL), 2, 19 , g,     A),
                                           (CB3R2R3e(x, q, :NORMAL), 3, 5.5, g,     A),
                                           (CB3R2R3c(x, q, :NORMAL), 3, 5.8, g,     A)]
                                           # (CB4R3R4( x, q, :NORMAL), 4, 0.16)]

        # exact values of integral
        exact = [5, exp(5) - exp(0), 5^2/2]

        # error should decrease at certain rate
        for Δt = logspace(0, -2.5, 5)

            # forward map
            f = flow(_g, _A, quad, method, TimeStepConstant(Δt))

            # initial conditions
            x₀ = Float64[1.0]
            q₀ = Float64[0.0, 0.0, 0.0]

            # call
            f(x₀, q₀, (0, 5), mon)

            # test 
            @test eltype(samples(mon)) == Coupled{Vector{Float64}, Vector{Float64}}

            # integrals
            @test norm(abs.(q₀ - exact)) / Δt^order < value

            # test allocations
            fun(f, x₀, q₀, span) = @allocated f(x₀, q₀, span)
            @test fun(f, x₀, q₀, (0.0, 100.0)) <= 32
        end
    end
end