@testset "test RAMStorage                        " begin
    store = RAMStorage{Float64, Tuple{Float64, Int}}()

    push!(store, 0.0, (0.0, 0))
    push!(store, 1.0, (1.0, 1))
    push!(store, 2.0, (2.0, 2))

    @test times(store)   == [0.0, 1.0, 2.0]
    @test samples(store) == [(0.0, 0), (1.0, 1), (2.0, 2)]

    reset!(store, 0)
    @test length(times(store))   == 0
    @test length(samples(store)) == 0
end

@testset "cubic interpolation                    " begin
    # data
    ts = 0.0:0.1:1

    # create storage with one sample
    sol1 = RAMStorage{Float64, Vector{Float64}}()
    sol2 = RAMStorage{Float64, Vector{Float64}}()
    sol3 = RAMStorage{Float64, Vector{Float64}}()

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

    	# check derivative (why do I want this????)
        @test norm(sol1(out, ti, Val{1}()) - [1]) < 1e-14
        @test norm(sol2(out, ti, Val{1}()) - [2*ti]) < 1e-14
        @test norm(sol3(out, ti, Val{1}()) - [3*ti*ti]) < 1e-14
    end

    # check out of bound range
    @test_throws ErrorException sol1(out, -1.0)
    @test_throws ErrorException sol1(out,  1.1)
end