@testset "stage cache tests                      " begin 
    NS = 1
    for c in (RAMStageCache(NS, Vector{Float64}), RAMStageCache(NS, [0.0]))
        push!(c, 0.0, 0.1, ([0.0],))
        push!(c, 0.1, 0.1, ([1.0],))
        push!(c, 0.2, 0.1, ([2.0],))

        @test c.ts  == [0.0, 0.1, 0.2]
        @test c.Δts == [0.1, 0.1, 0.1]
        @test c.xs  == [([0.0],), ([1.0],), ([2.0],)]
        
        reset!(c)
        @test c.ts  == []
        @test c.Δts == []
        @test c.xs  == []

        d = similar(c)
        @test d.ts  == []
        @test d.Δts == []
        @test d.xs  == []
    end
end