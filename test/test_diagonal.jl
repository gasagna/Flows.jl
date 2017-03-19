using Base.Test
using IMEXRK

@testset "verify interface" begin
    
    g(t, x, ẋ) = (ẋ .= 0.5*x; ẋ)
    A = Diagonal([0.5])
    f = DiagonalIMEXSystem(g, A)

    # 
    @test stiff(f) == A
    @test nonstiff(f) == g

    # 
    @test A_mul_B!(A, [1.0], [0.0]) == [0.5]

    # 
    @test ImcA!(stiff(f), 0.1, [1.0], [0.0]) == [1/0.95]
end