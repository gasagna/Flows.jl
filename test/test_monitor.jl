import LinearAlgebra: Diagonal, norm

@testset "test monitor type                      " begin
    m = Monitor(0, string)
    push!(m, 0.0, 0)
    @test times(m)   == [0.0]
    @test samples(m) == ["0"]
    @test eltype(samples(m)) == String
end

@testset "test monitor content                   " begin
	# store every push
    m = Monitor([1.0], x->x[1]^2)

    push!(m, 0.0, [0.0])
    push!(m, 0.1, [1.0])
    push!(m, 0.2, [2.0])
    push!(m, 0.3, [3.0])

    @test times(m)   == [0.0, 0.1, 0.2, 0.3]
    @test samples(m) == [0.0, 1.0, 4.0, 9.0]

    # skip one sample
    m = Monitor([1.0], x->x[1]^2; oneevery=2)

    push!(m, 0.0, [0.0])
    push!(m, 0.1, [1.0])
    push!(m, 0.2, [2.0])
    push!(m, 0.3, [3.0])

    @test times(m)   == [0.0, 0.2]
    @test samples(m) == [0.0, 4.0]
end

@testset "storeonebutlast                        " begin
    # integral of t in dt
    g(t, x, dxdt) = (dxdt[1] = t; dxdt)
    A = Diagonal([0.0])

    # integration scheme
    scheme = CB3R2R3e(Float64[0.0], :NORMAL)

    # monitors
    m = StoreOneButLast(zeros(1), identity)

    # forward map
    ϕ = flow(g, A, scheme, TimeStepConstant(0.1))

    # initial condition
    x₀ = [0.0]

    # test end point is calculated correctly
    ϕ(x₀, (0, 1), m)

    @test m.t == 0.9
    @test m.x == [0.5*(0.9)^2]
end

@testset "allocation                             " begin
    # integral of t in dt
    g(t, x, ẋ) = (ẋ[1] = t; ẋ)
    A = Diagonal([0.0])

    # integration scheme
    scheme = CB3R2R3e(Float64[0.0], :NORMAL)

    # monitors
    m = Monitor([1.0], x->x[1]^2; sizehint=10000)

    # forward map
    ϕ = flow(g, A, scheme, TimeStepConstant(0.01))

    # initial condition
    x₀ = [0.0]

    # test end point is calculated correctly
    ϕ(x₀, (0, 1.005), m)

    @test times(m)[end  ] == 1.005
    @test times(m)[end-1] == 1.000
    @test times(m)[end-2] == 0.990

    # warm up
    fun(ϕ, x₀, span, m) = @allocated ϕ(x₀, span, m)

    # does not allocate because we do not grow the arrays in the monitor
    @test fun(ϕ, x₀, (0, 1), m) == 0

    # try resetting and see we still do not allocate
    reset!(m)
    @test (@allocated reset!(m)) == 0

    @test fun(ϕ, x₀, (0, 1), m) == 0
end

@testset "reset monitors                         " begin
    m = Monitor([1.0], x->x[1]^2)

    push!(m, 0.0, [0.0])
    push!(m, 0.1, [1.0])
    push!(m, 0.2, [2.0])
    push!(m, 0.3, [3.0])

    @test times(m)   == [0.0, 0.1, 0.2, 0.3]
    @test samples(m) == [0.0, 1.0, 4.0, 9.0]

    reset!(m)
    @test times(m)   == []
    @test samples(m) == []
end

@testset "savebetween                            " begin
    m = Monitor([0.0], copy; savebetween=(1, 2))
    push!(m, 0, [0.0])
    push!(m, 1, [0.0])
    push!(m, 2, [0.0])
    push!(m, 3, [0.0])
    @test length(times(m)) == 2
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
    	# check function value
        @test norm(sol1(out, ti, Val{0}()) - [ti]) < 1e-14
        @test norm(sol2(out, ti, Val{0}()) - [ti*ti]) < 1e-14
        @test norm(sol3(out, ti, Val{0}()) - [ti*ti*ti]) < 1e-14

    	# check derivative
        @test norm(sol1(out, ti, Val{1}()) - [1]) < 1e-14
        @test norm(sol2(out, ti, Val{1}()) - [2*ti]) < 1e-14
        @test norm(sol3(out, ti, Val{1}()) - [3*ti*ti]) < 1e-14
    end

    # check out of bound range
    @test_throws ErrorException sol1(out, -1.0)
    @test_throws ErrorException sol1(out,  1.1)
end