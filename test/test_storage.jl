using Base.Test
using IMEXRKCB

@testset "test RAMStorage                        " begin
    store = RAMStorage{Float64, Tuple{Float64, Int}}()

    push!(store, 0.0, (0.0, 0))
	push!(store, 1.0, (1.0, 1))
	push!(store, 2.0, (2.0, 2))

	@test times(store)   == [0.0, 1.0, 2.0]
	@test samples(store) == [(0.0, 0), (1.0, 1), (2.0, 2)]

	resize!(store, 0)
	@test length(times(store))   == 0
	@test length(samples(store)) == 0
end