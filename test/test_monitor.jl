@testset "test monitor content                   " begin
    m = Monitor([1.0], x->x[1]^2)

    push!(m, 0.0, [0.0])
    push!(m, 0.1, [1.0])
    push!(m, 0.2, [2.0])
    push!(m, 0.3, [3.0])

    @test m.ts == [0.0, 0.1, 0.2, 0.3]
    @test m.xs == [0.0, 1.0, 4.0, 9.0]
end

@testset "allocation                             " begin
    # integral of t in dt
    g(t, x, ẋ) = (ẋ[1] = t; ẋ)
    A = Diagonal([0.0])

    # integration scheme
    scheme = IMEXRKScheme(IMEXRK3R2R(IMEXRKCB3e, false), [0.0])

    # monitors
    m = Monitor([1.0], x->x[1]^2)

    # forward map
    ϕ = integrator(g, A, scheme, 0.1)
    
    # initial condition
    x₀ = [0.0]

    # warm up
    ϕ(x₀, (0, 1), m) 

    # does not allocate because we do not grow the arrays in the monitor
    @test (@allocated ϕ(x₀, (0, 1), m)) == 0

    # try resetting and see we still do not allocate
    reset!(m)
    @test (@allocated reset!(m)) == 0
    
    @test (@allocated ϕ(x₀, (0, 1), m)) == 0
end

@testset "reset monitors                         " begin
    m = Monitor([1.0], x->x[1]^2)
    
    push!(m, 0.0, [0.0])
    push!(m, 0.1, [1.0])
    push!(m, 0.2, [2.0])
    push!(m, 0.3, [3.0])

    @test m.ts == [0.0, 0.1, 0.2, 0.3]
    @test m.xs == [0.0, 1.0, 4.0, 9.0]

    reset!(m)
    @test m.ts == []
    @test m.xs == []
end

@testset "cubic interpolation                    " begin
    # data
    ts = 0.0:0.1:1

    # create storage with one sample
    sol = Monitor([0.0])

    # push a cubic function of time
    for t in ts
        push!(sol, t, [t*t*t])
    end

    # check interpolation
    out = [0.0]
    for ti in 0.0:0.01:1.0
        out = sol(out, ti)
        @test out ≈ [ti*ti*ti]
    end

    # check out of bound range
    @test_throws ErrorException sol(out, -1.0)
    @test_throws ErrorException sol(out,  1.1)
end