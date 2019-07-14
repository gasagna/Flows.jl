@testset "couple                                 " begin
    # couplecopy
    a = couplecopy(2, [1.0])
    @test a[1] == [1.0]
    @test a[2] == [1.0]
    @test length(a) == 2
end