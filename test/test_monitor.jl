using Base.Test
using IMEXRKCB

@testset "monitor      " begin

    f1(x) = x[1]
    m = Monitor([1.0], f1)
    
    push!(m, [0.0], 0.0)
    push!(m, [1.0], 0.1)
    push!(m, [2.0], 0.2)
    push!(m, [3.0], 0.3)

    @test m.times   == [0.0, 0.1, 0.2, 0.3]
    @test m.samples == [0.0, 1.0, 2.0, 3.0]
end