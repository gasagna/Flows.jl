@testset "integrator                             " begin
    # make system
    g(t, x, ẋ) = (ẋ .= -0.5.*x; ẋ)
    A = Diagonal([-0.5])

    # define state
    x = [1.0]

    # integration scheme
    for impl in [IMEXRK3R2R(IMEXRKCB3e, false),
                 IMEXRK3R2R(IMEXRKCB3c, false),
                 IMEXRK4R3R(IMEXRKCB4,  false)]

        # define scheme
        scheme = IMEXRKScheme(impl, x)                 

        # integrator
        I = integrator(g, A, scheme, 0.01123)

        # map generator
        gen = fwdmapgen(I)

        # example
        for T = [0.1, 0.2, 0.5]
            f = gen(T)

            # test value
            @test f([1.0]) ≈ [exp(-T)]

            # test no allocations
            x₀ = Float64[1.0]
            @test (@allocated f(x₀)) == 0
        end
    end
end