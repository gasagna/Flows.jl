using Base.Test
using Flows

# example time stepping from hook
struct TestFromHook <: AbstractTimeStepFromHook end
(::TestFromHook)(g, A, z) = 0.1*sqrt(z[1])

# system
g(t, x, ẋ) = (ẋ .= -1; ẋ)

@testset "from hook                              " begin
    # integration scheme
    scheme = RK4(zeros(1), :NL)

    # forward map
    ϕ = flow(g, scheme, TestFromHook())

    # monitor
    mon = Monitor(zeros(1), x->x[1])

    # FORWARD INTEGRATION
    # run
    ϕ([1.0], (0, 1), reset!(mon))

    # test we hit the end point
    @test times(mon)[1]   == 0
    @test times(mon)[end] == 1

    # and that the solution is OK
    @test abs(samples(mon)[end]) < 2e-16

    # backward integration fails
    @test_throws ArgumentError ϕ([1.0], (1, 0), reset!(mon))
end

@testset "constant time step                     " begin
    # integration scheme
    scheme = RK4(zeros(1), :NL)

    # forward map
    ϕ = flow(g, scheme, TimeStepConstant(0.5))

    # monitor
    mon = Monitor(zeros(1), x->x[1])

    # FORWARD INTEGRATION
    # run
    ϕ([1.0], (0, 1), reset!(mon))

    # test we hit the end point
    @test times(mon)[1] == 0.0
    @test times(mon)[2] == 0.5
    @test times(mon)[3] == 1.0

    # and that the solution is OK
    @test abs(samples(mon)[end]) < 2e-16

    # backward integration fails
    @test_throws ArgumentError ϕ([1.0], (1, 0), reset!(mon))
end
