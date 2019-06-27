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

@testset "_interp_range                          " begin
    # data with period 1.0. Note t must always be in bounds
    #     1    2    3    4    5    6    7    8    9    10
    ts = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    # periodic case: note the difference between t = 0.0 and t = 1.0
    @test Flows._interp_indices( 0.00, ts, Val(4), true) == (10,  1, 2, 3)
    @test Flows._interp_indices( 0.01, ts, Val(4), true) == (10,  1, 2, 3)
    @test Flows._interp_indices( 0.11, ts, Val(4), true) == ( 1,  2, 3, 4)
    @test Flows._interp_indices( 0.51, ts, Val(4), true) == ( 5,  6, 7, 8)
    @test Flows._interp_indices( 0.91, ts, Val(4), true) == ( 9, 10, 1, 2)
    @test Flows._interp_indices( 1.00, ts, Val(4), true) == ( 9, 10, 1, 2)

    # non periodic case
    @test Flows._interp_indices(0.01, ts, Val(4), false) == (1, 2, 3, 4)
    @test Flows._interp_indices(0.11, ts, Val(4), false) == (1, 2, 3, 4)
    @test Flows._interp_indices(0.51, ts, Val(4), false) == (5, 6, 7, 8)
    @test Flows._interp_indices(0.51, ts, Val(6), false) == (4, 5, 6, 7, 8, 9)
    @test Flows._interp_indices(0.69, ts, Val(4), false) == (6, 7, 8,  9)
    @test Flows._interp_indices(0.71, ts, Val(4), false) == (7, 8, 9, 10)
    @test Flows._interp_indices(0.81, ts, Val(4), false) == (7, 8, 9, 10)
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

@testset "interpolation - non-periodic           " begin
    # data
    ts = 0.0:0.1:1

    # use these degrees for the interpolating polynomials
    degrees = [3, 5, 7, 9]

    # create storages of increasing order
    sols = ntuple(j->RAMStorage(Vector{Float64};
                                degree=degrees[j]), length(degrees))

    # number of interpolation points be smaller than length of data
    for (i, sol) in enumerate(sols)
        @test_throws ArgumentError sol(zeros(1), 0.1, Val(0))
    end

    # push a function of time that the interpolator should follow exactly
    for t in ts
        for (i, sol) in enumerate(sols)
            push!(sol, t, [t^degrees[i]])
        end
    end

    # check interpolation on finer grid
    out = [0.0]
    for (i, sol) in enumerate(sols)
        @test_throws ArgumentError sol(out, -0.1, Val(0))
        @test_throws ArgumentError sol(out,  1.1, Val(0))
        for ti in 0.0:0.01:1.0
            # check function value
            @test norm(sol(out, ti, Val(0)) - [ti^degrees[i]]) < 1e-14
        end
        # check derivative (why do I want this???? when you need interpolation of )
    end
end

@testset "_wrap_around_point                     " begin
    @test  Flows._wrap_around_point((  1,   2,   3, 4)) == 5
    @test  Flows._wrap_around_point((100,   1,   2, 3)) == 1
    @test  Flows._wrap_around_point(( 99, 100,   1, 2)) == 2
    @test  Flows._wrap_around_point(( 98,  99, 100, 1)) == 3
end

@testset "interpolation - periodic               " begin
    # data
    ts = range(0, stop=2.0, length=101)[1:100]

    # use these degrees for the interpolating polynomials
    degrees = [3, 5, 7, 9]

    # create storages of increasing order
    sols = ntuple(j->RAMStorage(Vector{Float64};
                                degree=degrees[j],
                                period=2.0), length(degrees))

    # tolerances depend on order
    tols = [1e-6, 1e-9, 1e-12, 1e-15]

    # push a function of time that the interpolator should follow exactly
    for t in ts
        for (i, sol) in enumerate(sols)
            push!(sol, t, Float64[cos(π*t)])
        end
    end

    #
    fun(sol, out) = @allocated sol(out, 0.5, Val(0))

    # check interpolation on another grid
    out = [0.0]
    for i in 1:length(degrees)
        # test allocations
        @test fun(sols[i], out) == 0

        # test time out of bounds 
        @test_throws ArgumentError sols[i](out, -0.1, Val(0))
        @test_throws ArgumentError sols[i](out,  2.1, Val(0))
        for ti in range(0.0, stop=2.0, length=17)
            # check function value
            @test norm(sols[i](out, ti, Val(0)) - Float64[cos(π*ti)]) < tols[i]
            # check derivative (why do I want this???? )
        end
    end
end