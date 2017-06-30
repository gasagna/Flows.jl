@testset "diffusion                              " begin
    g(t, x, ẋ) = (ẋ[1] = 1)
    A = Diagonal([-5])

    # integration scheme
    for scheme in [IMEXRK3R2R(IMEXRKCB3e, false, [0.0]),
                   IMEXRK3R2R(IMEXRKCB3c, false, [0.0]),
                   IMEXRK4R3R(IMEXRKCB4,  false, [0.0])]

        # monitor
        m = Monitor([0.0], x->x[1])

        # forward map
        ϕ = integrator(g, A, scheme, 0.12345, m)
        
        # do job
        ϕ([0.0], 20)

        # 
        @test m.samples[end] ≈ 1/5
    end
end