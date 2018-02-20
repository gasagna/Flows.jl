#=
    Comparing the numerical solution with a case for which there is a
    # known analytic solution, discussed in the tutorial at 
    http://cs.stanford.edu/~ambrad/adjoint_tutorial.pdf, section 2.1.3.
=#

# define adjoint equation
struct AdjointEq{M<:Monitor}
    mon::M       # storage for forward solution
      b::Float64 # coefficient
end

# integrate this fully explicitly
(eq::AdjointEq)(t::Real, λ, dλdt) = (dλdt .= 1 .- eq.b.*λ)

# define quadrature function for gradient integration
struct QuadFun{X, M<:Monitor}
    mon::M # storage for forward solution
    tmp::X # temporary to calculate state
end

# explicit integration of the quadrature rule (note negative sign)
(qfun::QuadFun)(t::Real, λ, dqdt) = dqdt .= .-(.- λ .* qfun.mon(qfun.tmp, t))

@testset "adjoint                                " begin
    # constant parameters
    const a, b, T = 1.0, 2.0, 1.0

    # define forward system, fully explicit
    g_forw(t, x, dxdt) = dxdt .= b.*x

    # initial condition
    x0   = Float64[a]
    λ0   = Float64[0.0]
    grad = Float64[0.0]

    # define storage for forward solution and monitor for backward
    forw_sol = Monitor(x0, copy)
    back_sol = Monitor(λ0, copy)

    # instantiate objects
    g_back  = AdjointEq(forw_sol, b)
    quadfun = QuadFun(forw_sol, similar(grad))

    # define forward and backward integrators
    forw_method = IMEXMethod(:CB4_4R3R, x0)
    back_method = IMEXMethod(:CB4_4R3R, λ0, grad)
    G_forw = integrator(g_forw, forw_method, 0.01)
    G_back = integrator(g_back, nothing, quadfun, back_method, 0.01)

    # solve forward problem
    G_forw(x0, (0, T), forw_sol)

    # check forward solution
    exact = a*exp.(b*forw_sol.ts)
    simul = [x[1] for x in forw_sol.xs]
    @test maximum(abs.(exact - simul)) < 1e-8

    # solve backward problem
    G_back(λ0, grad, (T, 0), back_sol)

    # check backward problem solution
    exact = (1 - exp.(b*(T - back_sol.ts)))./b
    simul = [x[1] for x in back_sol.xs]
    @test maximum(abs.(exact - simul)) < 1e-8

    # check gradient
    simul = grad[1]
    exact = a/b*T*exp(b*T) - a/b^2*(exp(b*T) - 1)
    @test abs(exact - simul) < 1e-8
end