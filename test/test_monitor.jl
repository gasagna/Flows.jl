using Base.Test
using IMEXRKCB

@testset "one monitor  " begin
    m = Monitor([1.0], x->x[1])
    
    push!(m, [0.0], 0.0)
    push!(m, [1.0], 0.1)
    push!(m, [2.0], 0.2)
    push!(m, [3.0], 0.3)

    @test m.times   == [0.0, 0.1, 0.2, 0.3]
    @test m.samples == [0.0, 1.0, 2.0, 3.0]
end


@testset "two monitors " begin
    m1 = Monitor([1.0], x->x[1]  )
    m2 = Monitor([1.0], x->x[1]^2)
    
    push!((m1, m2), [0.0], 0.0)
    push!((m1, m2), [1.0], 0.1)
    push!((m1, m2), [2.0], 0.2)
    push!((m1, m2), [3.0], 0.3)

    @test m1.times   == [0.0, 0.1, 0.2, 0.3]
    @test m2.times   == [0.0, 0.1, 0.2, 0.3]
    @test m1.samples == [0.0, 1.0, 2.0, 3.0]
    @test m2.samples == [0.0, 1.0, 4.0, 9.0]
end