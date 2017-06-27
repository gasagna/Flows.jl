@testset "verify interface                       " begin
    A = Diagonal([0.5])
    @test ImcA!(A, 0.1, [1.0], [0.0]) == [1/0.95]
end