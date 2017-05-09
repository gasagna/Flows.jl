using Base.Test
using IMEXRKCB


@testset "diffusion    " begin
    g(t, x, ẋ) = (ẋ[1] = 1)
    A = Diagonal([-5])

    # integration scheme
    for scheme in [IMEXRK3R2R(IMEXRKCB3e, [0.0], false),
                   IMEXRK3R2R(IMEXRKCB3c, [0.0], false),
                   IMEXRK4R3R(IMEXRKCB4,  [0.0], false)]

        # forward map
        ϕ = forwmap!(g, A, 20, 0.12345, scheme)

        # monitor
        m = Monitor([0.0], x->x[1])
        
        # do job
        ϕ([0.0], m)

        # 
        @test m.samples[end] ≈ 1/5
    end
end