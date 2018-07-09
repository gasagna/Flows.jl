@testset "diffusion                              " begin
    g(t, x, ẋ) = (ẋ[1] = 1)
    A = Diagonal([-5])

    # integration scheme
    for method in [Scheme(:CB2_3R2R,  Float64[1.0]),
                   Scheme(:CB3e_3R2R, Float64[1.0]),
                   Scheme(:CB3c_3R2R, Float64[1.0]),
                   Scheme(:CB4_4R3R,  Float64[1.0])]

        # monitor
        m = Monitor([0.0], x->x[1])

        # forward map
        ϕ = integrator(g, A, method, TimeStepConstant(0.12345))
        
        # do job
        ϕ([0.0], (0, 20), m)

        # test value
        @test samples(m)[end] ≈ 1/5

        # ensure we saved the last step in the monitor
        @test times(m)[end] == 20
    end
end