@testset "losslessrange                          " begin
    @testset "lossy                          " begin
        rng = IMEXRKCB.LossLessRange(0, 1, 0.4)
        expected = [0, 0.4, 0.8, 1.0]
        @test collect(rng) == expected
        for (i, el) in enumerate(rng)
            @test el == expected[i]
        end
        @test length(rng) == 4
        @test rng[1] === 0.0
        @test rng[2] === 0.4
        @test rng[3] === 0.8
        @test rng[4] === 1.0
        @test_throws BoundsError rng[0]
        @test_throws BoundsError rng[5]
    end
    @testset "lossless                           " begin
        rng = IMEXRKCB.LossLessRange(0, 0.8, 0.4)
        expected = [0, 0.4, 0.8]
        @test collect(rng) == expected
        for (i, el) in enumerate(rng)
            @test el == expected[i]
        end
        @test length(rng) == 3
        @test rng[1] === 0.0
        @test rng[2] === 0.4
        @test rng[3] === 0.8
        @test_throws BoundsError rng[0]
        @test_throws BoundsError rng[4]
    end
    @testset "backwards lossy                    " begin
        rng = IMEXRKCB.LossLessRange(1, 0, -0.4)
        expected = [1.0, 0.6, 0.2, 0.0]
        @test collect(rng) == expected
        for (i, el) in enumerate(rng)
            @test el == expected[i]
        end
        @test length(rng) == 4
        @test rng[1] === 1.0
        @test rng[2] === 0.6
        @test rng[3] === 0.2
        @test rng[4] === 0.0
        @test_throws BoundsError rng[0]
        @test_throws BoundsError rng[5]
    end
    @testset "backwards lossless                    " begin
        rng = IMEXRKCB.LossLessRange(0.8, 0, -0.4)
        expected = [0.8, 0.4, 0.0]
        @test collect(rng) == expected
        for (i, el) in enumerate(rng)
            @test el == expected[i]
        end
        @test length(rng) == 3
        @test rng[1] === 0.8
        @test rng[2] === 0.4
        @test rng[3] === 0.0
        @test_throws BoundsError rng[0]
        @test_throws BoundsError rng[4]
    end
end
