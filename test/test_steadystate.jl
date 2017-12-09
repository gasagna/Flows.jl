@testset "diffusion                              " begin
    g(t, x, ẋ) = (ẋ[1] = 1)
    A = Diagonal([-5])

    # integration scheme
    for impl in [IMEXRK3R2R(IMEXRKCB2,  false),
                 IMEXRK3R2R(IMEXRKCB3e, false),
                 IMEXRK3R2R(IMEXRKCB3c, false),
                 IMEXRK4R3R(IMEXRKCB4,  false)]

        # define scheme                
        scheme = IMEXRKScheme(impl, [1.0])

        # monitor
        m = Monitor([0.0], x->x[1])

        # forward map
        ϕ = integrator(g, A, scheme, 0.12345)
        
        # do job
        ϕ([0.0], (0, 20), m)

        # test value
        @test m.xs[end] ≈ 1/5

        # ensure we saved the last step in the monitor
        @test m.ts[end] == 20
    end
end