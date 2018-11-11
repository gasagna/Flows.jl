import LinearAlgebra: Diagonal, norm

@testset "coupled state broadcast                " begin
    x = [1,   2,   3]
    q = [3.0, 4.0, 5.0]

    # define two coupled states
    z = Flows.couple(x, q)
    y = Flows.couple(copy(x).+1, copy(q).+1)

    # define some operation
    fun!(z, y) = (@all z .= 2.0.*z .+ 1.0.*y .- 1; z)

    # apply 
    fun!(z, y)

    # value
    @test z[1] == [3, 6, 9]
    @test z[2] == [9.0, 12.0, 15.0]

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
    @inline function quad(t, x, dxdt, q, dqdt)
        dqdt[1] = 1
        dqdt[2] = x[1]
        dqdt[3] = t
        return dqdt
    end

    # state and quadrature
    x, q = Float64[0.0], Float64[0.0, 0.0, 0.0]

    # monitor first quadrature equation
    mon = Monitor(couple(x, q), copy)

    # integration method
    for (method, order, value, _g, _A) in [(RK4(     couple(x, q), :NORMAL), 4, 6.2, gfull, nothing),
                                           (CB3R2R2( couple(x, q), :NORMAL), 2, 19 , g,     A),
                                           (CB3R2R3e(couple(x, q), :NORMAL), 3, 5.5, g,     A),
                                           (CB3R2R3c(couple(x, q), :NORMAL), 3, 5.8, g,     A)]
                                           # (CB4R3R4( x, q, :NORMAL), 4, 0.16)]

        # exact values of integral
        exact = [5, exp(5) - exp(0), 5^2/2]

        # error should decrease at certain rate
        for Δt = 10 .^ range(0, stop=-2.5, length=5)

            # forward map
            f = flow(couple(_g, quad), couple(_A, nothing), method, TimeStepConstant(Δt))

            # initial conditions
            x₀ = Float64[1.0]
            q₀ = Float64[0.0, 0.0, 0.0]

            # call
            f(couple(x₀, q₀), (0, 5), mon)

            # test 
            @test eltype(samples(mon)) == Coupled{2, Tuple{Vector{Float64}, Vector{Float64}}}

            # integrals
            @test norm(abs.(q₀ - exact)) / Δt^order < value

            # test allocations
            fun(f, xq, span) = @allocated f(xq, span)
            @test fun(f, couple(x₀, q₀), (0.0, 1000.0)) <= 32
        end
    end
end