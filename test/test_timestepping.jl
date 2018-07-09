using Base.Test
using Flows

# example time stepping from hook
struct TestFromHook <: AbstractTimeStepFromHook end
(::TestFromHook)(g, A, z) = 0.1*sqrt(z[1])

@testset "from hook                              " begin
    # make system
    g(t, x, ẋ) = (ẋ .= -1; ẋ)

    # integration scheme
    scheme = Scheme(:RK4, zeros(1))

    # forward map
    ϕ = integrator(g, nothing, scheme, TestFromHook())

    # monitor
    mon = Monitor(zeros(1), x->x[1])

    # FORWARD INTEGRATION
    # run
    ϕ([1.0], (0, 1), mon)

    # test we hit the end point
    @test times(mon)[1]   == 0
    @test times(mon)[end] == 1

    # and that the solution is OK
    @test abs(samples(mon)[end]) < 2e-16

    # BACKWARD INTEGRATION
    # run
    ϕ([1.0], (1, 0), reset!(mon))

    # test we hit the end point
    @test times(mon)[1]   == 1
    @test times(mon)[end] == 0

    # and that the solution is OK
    @test abs(samples(mon)[end] - 2.0) < 5e-16
end
