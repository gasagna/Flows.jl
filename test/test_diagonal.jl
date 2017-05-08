@testset "verify interface" begin
    
    g(t, x, ẋ) = (ẋ .= 0.5*x; ẋ)
    A = Diagonal([0.5])

    @test A_mul_B!(A, [1.0], [0.0]) == [0.5]
    @test ImcA!(A, 0.1, [1.0], [0.0]) == [1/0.95]
end