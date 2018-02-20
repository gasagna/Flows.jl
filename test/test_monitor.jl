using Base.Test
using IMEXRKCB

@testset "test monitor type                      " begin
    m = Monitor([1.0], x->"1")
    push!(m, 0.0, [0.0])
    @test m.xs == ["1"]

    @test eltype(m.xs) == String
end

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

    # integration method
    method = IMEXMethod(:CB3e_3R2R, Float64[0.0])

    # monitors
    m1 = Monitor([1.0], x->x[1]^2; sizehint=10000)
    m2 = Monitor([1.0], x->x[1]^3; sizehint=10000)

    # forward map
    ϕ = integrator(g, A, method, 0.01)
    
    # initial condition
    x₀ = [0.0]

    # test end point is calculated correctly
    ϕ(x₀, (0, 1.005), m1, m2)

    @test m1.ts[end  ] == 1.005
    @test m1.ts[end-1] == 1.000
    @test m1.ts[end-2] == 0.990

    # warm up
    fun(ϕ, x₀, span, m1, m2) = @allocated ϕ(x₀, span, m1, m2)

    # does not allocate because we do not grow the arrays in the monitor
    @test fun(ϕ, x₀, (0, 1), m1, m2) == 0

    # try resetting and see we still do not allocate
    for m = [m1, m2]
        reset!(m)
        @test (@allocated reset!(m)) == 0
    end
    
    @test fun(ϕ, x₀, (0, 1), m1, m2) == 0
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
    sol1 = Monitor([0.0])
    sol2 = Monitor([0.0])
    sol3 = Monitor([0.0])

    # push a cubic function of time
    for t in ts
        push!(sol1, t, [t])
        push!(sol2, t, [t*t])
        push!(sol3, t, [t*t*t])
    end

    # check interpolation on finer grid
    out = [0.0]
    for ti in 0.0:0.01:1.0
        @test sol1(out, ti) ≈ [ti]
        @test sol2(out, ti) ≈ [ti*ti]
        @test sol3(out, ti) ≈ [ti*ti*ti]
    end

    # check out of bound range
    @test_throws ErrorException sol1(out, -1.0)
    @test_throws ErrorException sol1(out,  1.1)
end