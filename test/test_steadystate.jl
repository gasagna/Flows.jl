using Base.Test
using IMEXRKCB

@testset "diffusion                              " begin
    g(t, x, ẋ) = (ẋ[1] = 1)
    A = Diagonal([-5])

    # integration scheme
    for scheme in [IMEXRK3R2R(IMEXRKCB3e, false, [0.0]),
                   IMEXRK3R2R(IMEXRKCB3c, false, [0.0]),
                   IMEXRK4R3R(IMEXRKCB4,  false, [0.0])]

        # monitor
        m = Monitor((x->x[1], ), [0.0])

        # forward map
        ϕ = integrator(g, A, scheme, 0.12345)
        
        # do job
        ϕ([0.0], 20, m)

        # 
        @test m.samples[1][end] ≈ 1/5
    end
end