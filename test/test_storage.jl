using Base.Test
using IMEXRKCB

# @testset "Test findsortidx                       " begin
#     x = [1.0, 2.0, 3.0, 4.0, 5.0]
#     @test IMEXRKCB.findsortidx(x, 3.2) == 3
#     @test IMEXRKCB.findsortidx(x, 1.0) == 1
#     @test IMEXRKCB.findsortidx(x, 2.0) == 2
#     @test IMEXRKCB.findsortidx(x, 5.0) == 5
#     @test IMEXRKCB.findsortidx(x, 1.0) == 1

#     x = [1.0, 2.0, 3.0, 4.0]
#     @test IMEXRKCB.findsortidx(x, 3.2) == 3
#     @test IMEXRKCB.findsortidx(x, 1.0) == 1
#     @test IMEXRKCB.findsortidx(x, 2.0) == 2
#     @test IMEXRKCB.findsortidx(x, 4.0) == 4
#     @test IMEXRKCB.findsortidx(x, 1.0) == 1

#     # return 0
#     x = 1:10
#     @test IMEXRKCB.findsortidx(x,  -1) == 0
#     @test IMEXRKCB.findsortidx(x,   0) == 0
#     @test IMEXRKCB.findsortidx(x,  11) == 0
#     @test IMEXRKCB.findsortidx(x,  5, 6, 3) == 0
#     @test IMEXRKCB.findsortidx(x,  5, 7, 9) == 0
#     @test IMEXRKCB.findsortidx(x,  9, 5, 6) == 0
# end

# @testset "Test functionality                     " begin
#     # data
#     N = 2
#     ts = 0.0:0.1:1

#     # create storage with one sample
#     sol = IMEXRKCB.SolutionStorage(zeros(N))

#     # push a cubic function of time
#     for t in ts
#         push!(sol, [t*t*t*j for j = 1:N], t)
#     end

#     # check interpolation
#     out = zeros(N)
#     for ti in 0.0:0.01:1.0
#         out = sol(out, ti)
#         @test out ≈ [ti*ti*ti, 2*ti*ti*ti]
#     end

#     # check out of bound range. Give a 0 idxcurr here, does not matter
#     @test_throws ErrorException sol(out, -1.0)
#     @test_throws ErrorException sol(out,  1.1)

#     # storage must be writable
#     @test_throws ErrorException push!(setrmode!(sol), zeros(N), 0)
#     @test iswmode(setrmode!(sol)) == false
#     @test iswmode(setwmode!(sol)) == true
#     @test isrmode(setrmode!(sol)) == true
#     @test isrmode(setwmode!(sol)) == false

#     # test resetting
#     reset!(sol)
#     @test length(sol.ts) == 0
#     @test length(sol.xs) == 0
#     @test sol.idx == 0
# end

# @testset "Test storage no allocations            " begin
#     # define system
#     g(t, x, ẋ) = (ẋ .= .-0.5.*x; ẋ)
#     A = Diagonal([-0.5])

#     # define integrator
#     scheme = IMEXRKScheme(IMEXRK3R2R(IMEXRKCB2, false), [1.0])
#     ϕ = integrator(g, A, scheme, 0.1)

#     # define initial condition and storage
#     x0 = [1.0]
#     storage = setwmode!(SolutionStorage(x0))

#     # warm up
#     ϕ(x0, 1, storage)

#     # test for no allocations
#     @test (@allocated ϕ(x0, 1, storage)) == 0
# end

# @testset "Test storage content                   " begin
#     # define system
#     g(t, x, ẋ) = (ẋ .= .-0.5.*x; ẋ)
#     A = Diagonal([-0.5])

#     # define integrator
#     scheme = IMEXRKScheme(IMEXRK4R3R(IMEXRKCB4, false), [1.0])
#     ϕ = integrator(g, A, scheme, 0.01)

#     # define initial condition and storage
#     x0 = [1.0]
#     storage = setwmode!(SolutionStorage(x0, true))

#     # warm up
#     ϕ(x0, 1, storage)

#     # test time
#     @test storage.ts ≈ collect(0:0.01:1)

#     # extract solution
#     xs = [x[1] for x in storage.xs]
#     @test maximum(abs.(xs - exp.(-storage.ts))) < 1e-8
# end

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

# explicit integration of the quadrature rule (note negative sign because backwards)
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
    forw_scheme = IMEXRKScheme(IMEXRK4R3R(IMEXRKCB4, false), x0)
    back_scheme = IMEXRKScheme(IMEXRK4R3R(IMEXRKCB4, false), λ0, grad)
    G_forw = integrator(g_forw, forw_scheme, 0.01)
    G_back = integrator(g_back, nothing, quadfun, back_scheme, 0.01)

    # solve forward problem
    G_forw(x0, (0, T), forw_sol)

    # check forward solution
    exact = a*exp.(b*forw_sol.ts)
    simul = [x[1] for x in forw_sol.xs]
    @test maximum(abs.(exact - simul)) < 1e-8

    # solve backward problem
    println(grad)
    G_back(λ0, grad, (T, 0), back_sol)
    println(grad)

    # check backward problem solution
    exact = (1 - exp.(b*(T - back_sol.ts)))./b
    simul = [x[1] for x in back_sol.xs]
    @test maximum(abs.(exact - simul)) < 1e-8

    # check gradient
    simul = grad[1]
    exact = a/b*T*exp(b*T) - a/b^2*(exp(b*T) - 1)
    @test abs(exact - simul) < 1e-8
end