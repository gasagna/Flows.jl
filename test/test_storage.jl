@testset "test RAMStorage constructor            " begin
    el = (0.0, 0)
    store = RAMStorage(el)
    store = RAMStorage(typeof(el))

    push!(store, 0.0, (0.0, 0))
    push!(store, 1.0, (1.0, 1))
    push!(store, 2.0, (2.0, 2))

    @test times(store)   == [0.0, 1.0, 2.0]
    @test samples(store) == [(0.0, 0), (1.0, 1), (2.0, 2)]

    reset!(store, 0)
    @test length(times(store))   == 0
    @test length(samples(store)) == 0

    # other argument values
    store = RAMStorage(el; ttype=Float32)
    @test eltype(times(store)) == Float32
    @test degree(store) == 3

    store = RAMStorage(el; ttype=Float16, degree=7, period=5.4578)
    @test eltype(times(store)) == Float16
    @test degree(store) == 7
end

@testset "_normalise                             " begin
    # periodic case
    @test Flows._normalise( 0.0, 1.0, true) ≈ 0.0
    @test Flows._normalise( 0.9, 1.0, true) ≈ 0.9
    @test Flows._normalise( 1.0, 1.0, true) ≈ 0.0
    @test Flows._normalise(-0.1, 1.0, true) ≈ 0.9
    @test Flows._normalise(-1.1, 1.0, true) ≈ 0.9
    @test Flows._normalise(-2.1, 1.0, true) ≈ 0.9
    @test Flows._normalise(+1.1, 1.0, true) ≈ 0.1
    @test Flows._normalise(+2.1, 1.0, true) ≈ 0.1

    # non-periodic
    @test Flows._normalise( 0.0, 1.0, false) == 0.0
    @test Flows._normalise( 1.0, 1.0, false) == 1.0
    @test Flows._normalise( 2.0, 1.0, false) == 2.0
end

@testset "_interp_range                          " begin
    # data with period 1.0. Note t must always be in bounds
    #     1    2    3    4    5    6    7    8    9    10
    ts = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    # number of interpolation points be smaller than length of data
    @test_throws ArgumentError Flows._interp_range(0.00, ts, Val(20), 1.0)

    # periodic case
    @test Flows._interp_range( 0.00, ts, Val(4), 1.0) == (10,  1, 2, 3)
    @test Flows._interp_range( 0.01, ts, Val(4), 1.0) == (10,  1, 2, 3)
    @test Flows._interp_range( 0.11, ts, Val(4), 1.0) == ( 1,  2, 3, 4)
    @test Flows._interp_range( 0.51, ts, Val(4), 1.0) == ( 5,  6, 7, 8)
    @test Flows._interp_range( 0.91, ts, Val(4), 1.0) == ( 9, 10, 1, 2)
    @test Flows._interp_range( 1.00, ts, Val(4), 1.0) == (10,  1, 2, 3)

    # non periodic case
    @test_throws ArgumentError Flows._interp_range(-0.1, ts, Val(4), 0.0)
    @test_throws ArgumentError Flows._interp_range(+1.0, ts, Val(4), 0.0)

    @test Flows._interp_range(0.01, ts, Val(4), 0.0) == (1, 2, 3, 4)
    @test Flows._interp_range(0.11, ts, Val(4), 0.0) == (1, 2, 3, 4)
    @test Flows._interp_range(0.51, ts, Val(4), 0.0) == (5, 6, 7, 8)
    @test Flows._interp_range(0.51, ts, Val(6), 0.0) == (4, 5, 6, 7, 8, 9)
    @test Flows._interp_range(0.69, ts, Val(4), 0.0) == (6, 7, 8,  9)
    @test Flows._interp_range(0.71, ts, Val(4), 0.0) == (7, 8, 9, 10)
    @test Flows._interp_range(0.81, ts, Val(4), 0.0) == (7, 8, 9, 10)
end

@testset "_lagr_weights                          " begin
    # trivial case
    ts = (0.0, 1.0, 2.0, 3.0, 4.0, 5.0)
    @test Flows._lagr_weights(0.0, ts, Val(0)) == (1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    @test Flows._lagr_weights(1.0, ts, Val(0)) == (0.0, 1.0, 0.0, 0.0, 0.0, 0.0)
    @test Flows._lagr_weights(4.0, ts, Val(0)) == (0.0, 0.0, 0.0, 0.0, 1.0, 0.0)
    @test Flows._lagr_weights(5.0, ts, Val(0)) == (0.0, 0.0, 0.0, 0.0, 0.0, 1.0)

    # linear case is a weighted average
    ts = (0.0, 1.0)
    @test Flows._lagr_weights(0.4, ts, Val(0)) == (0.6, 0.4)
end

@testset "cubic interpolation                    " begin
    # data
    ts = 0.0:0.1:1

    # use these degrees for the interpolating polynomials
    degrees = [3, 5, 7, 9]

    # create storages of increasing order
    sols = ntuple(j->RAMStorage(Vector{Float64}; degree=degrees[j]), length(degrees))

    # push a function of time that the interpolator should follow exactly
    for t in ts
        for (i, sol) in enumerate(sols)
            push!(sol, t, [t^degrees[i]])
        end
    end

    # check interpolation on finer grid
    out = [0.0]
    for ti in 0.0:0.01:1.0
        # check function value
        for (i, sol) in enumerate(sols)
            @test norm(sol(out, ti, Val(0)) - [ti^degrees[i]]) < 1e-14
        end
        # check derivative (why do I want this???? when you need interpolation of )
    end

    #     @test norm(sol1(out, ti, Val(1)) - [1]) < 1e-14
    #     @test norm(sol2(out, ti, Val(1)) - [2*ti]) < 1e-14
    #     @test norm(sol3(out, ti, Val(1)) - [3*ti*ti]) < 1e-14

    # check out of bound range
    # @test_throws ErrorException sol1(out, -1.0)
    # @test_throws ErrorException sol1(out,  1.1)
end