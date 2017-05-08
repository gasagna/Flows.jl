@testset "forwmap!" begin
    # make system
    g(t, x, ẋ) = (ẋ .= -0.5.*x; ẋ)
    A = Diagonal([-0.5])

    # vector
    x = [1.0]

    # integration scheme
    for scheme in [IMEXRK3R2R(IMEXRKCB3e, x, false),
                   IMEXRK3R2R(IMEXRKCB3c, x, false),
                   IMEXRK4R3R(IMEXRKCB4,  x, false)]

        # forward map
        ϕ = forwmap!(g, A, 1.0, 0.01123, scheme)

        # check relative error
        @test ((ϕ(x) - exp(-1))/exp(-1))[1] < 4.95e-8
        @test ((ϕ(x) - exp(-2))/exp(-2))[1] < 4.95e-8
        @test ((ϕ(x) - exp(-3))/exp(-3))[1] < 4.95e-8
        @test ((ϕ(x) - exp(-4))/exp(-4))[1] < 4.95e-8
        @test ((ϕ(x) - exp(-5))/exp(-5))[1] < 4.95e-8

        # try giving a different input
        @test_throws MethodError ϕ([1])
    end
end

@testset "time step" begin
    @test IMEXRKCB.next_Δt(0.0, 1.0, 0.1) == 0.1
    @test IMEXRKCB.next_Δt(0.0, 1.0, 1.1) == 1.0
end