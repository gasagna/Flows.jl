using Base.Test
using IMEXRK

@testset "linear system" begin


    # make system
    g(t, x, ẋ) = (ẋ .= -0.5.*x; ẋ)
    A = Diagonal([-0.5])
    f = DiagonalIMEXSystem(g, A)

    # vector
    x = [1.0]

    # integration scheme
    scheme = IMEXRK3R2R(IMEXRKCB3e, x)

    # forward map
    ϕ = forwmap!(f, 1.0, 0.01, scheme)

    # check relative error
    @test ((ϕ(x) - exp(-1))/exp(-1))[1] < 4e-8
    @test ((ϕ(x) - exp(-2))/exp(-2))[1] < 4e-8
    @test ((ϕ(x) - exp(-3))/exp(-3))[1] < 4e-8
    @test ((ϕ(x) - exp(-4))/exp(-4))[1] < 4e-8
    @test ((ϕ(x) - exp(-5))/exp(-5))[1] < 4e-8
end

@testset "time step" begin
    @test IMEXRK.next_Δt(0.0, 1.0, 0.1) == 0.1
    @test IMEXRK.next_Δt(0.0, 1.0, 1.1) == 1.0
end