@testset "basic                                  " begin
    x = [0, 0]
    y = [1, 2]
    @over_i x[i] = 2 * y[i]
    @test x == [2, 4]
end

@testset "operators                              " begin
    x = [0, 0]
    y = [1, 2]
    @over_i x[i] += 2 * y[i]
    @test x == [2, 4]

    @over_i x[i] *= 2 * y[i]
    @test x == [4, 16]

    @over_i x[i] -= 2 * y[i]
    @test x == [2, 12]

    @over_i x[i] /= 2
    @test x == [1, 6]
end

@testset "no allocations                         " begin
    x = rand(10000)
    y = rand(10000)
    f(x, y) = @over_i x[i] = 2 * y[i]
    f(x, y)
    @test (@allocated f(x, y)) == 0
    @test x == 2*y
end