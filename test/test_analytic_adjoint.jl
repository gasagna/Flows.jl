# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #

#=
    Test the whole framework by comparing the numerical solution with 
    a simple case for which there is a known analytic solution, discussed
    in the tutorial at http://cs.stanford.edu/~ambrad/adjoint_tutorial.pdf
    in section 2.1.3.
=#

# The system is a 1D linear system 
struct ExampleSystem 
    b::Float64 # jacobian
    f::Float64 # forcing
end

# linear problem
(eq::ExampleSystem)(t, x, dxdt) = (dxdt[1] = eq.b*x[1] + eq.f; dxdt)

# adjoint problem
(eq::ExampleSystem)(t, x, q, dqdt) = (dqdt[1] = eq.b*q[1] + eq.f; dqdt)

@testset "analytic adjoint check                 " begin

    # parameters
    a, b, T = 0.67, 0.8, 1.23

    # initial condition
    x0 = Float64[a]

    # explicit integrator
    method = RK4(x0)

    # flow map
    ψ = flow(ExampleSystem(b, 0), method, TimeStepConstant(1e-3))

    # get solution and fill the cache used later for the adjoint calculations
    mon = Monitor(x0, x->x[1])
    cache = RAMStageCache(4, x0)
    ψ(copy(x0), (0, T), cache)
    ψ(copy(x0), (0, T), mon)

    # should be an exponential
    for (t, s) in zip(times(mon), samples(mon))
        @test abs(a*exp(b*t) - s) < 1e-13
    end

    # now solve adjoint problem

    # initial condition
    q0 = [0.0]

    # explicit integrator
    method = RK4(q0, true)

    # flow map (the system is self adjoint)
    ψ = flow(ExampleSystem(b, -1), method, TimeStepFromCache())

    # get solution
    ψ(q0, cache, reset!(mon))

    # should be an exponential
    for (t, s) in zip(times(mon), samples(mon))
        @test abs( 1/b*(1-exp(b*(T-t))) - s ) < 1e-13
    end
end