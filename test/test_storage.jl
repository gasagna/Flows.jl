using Base.Test
using IMEXRKCB

@testset "Test findsortidx                       " begin
    x = [1.0, 2.0, 3.0, 4.0, 5.0]
    @test IMEXRKCB.findsortidx(x, 3.2) == 3
    @test IMEXRKCB.findsortidx(x, 1.0) == 1
    @test IMEXRKCB.findsortidx(x, 2.0) == 2
    @test IMEXRKCB.findsortidx(x, 5.0) == 5
    @test IMEXRKCB.findsortidx(x, 1.0) == 1

    x = [1.0, 2.0, 3.0, 4.0]
    @test IMEXRKCB.findsortidx(x, 3.2) == 3
    @test IMEXRKCB.findsortidx(x, 1.0) == 1
    @test IMEXRKCB.findsortidx(x, 2.0) == 2
    @test IMEXRKCB.findsortidx(x, 4.0) == 4
    @test IMEXRKCB.findsortidx(x, 1.0) == 1

    # return 0
    x = 1:10
    @test IMEXRKCB.findsortidx(x,  -1) == 0
    @test IMEXRKCB.findsortidx(x,   0) == 0
    @test IMEXRKCB.findsortidx(x,  11) == 0
    @test IMEXRKCB.findsortidx(x,  5, 6, 3) == 0
    @test IMEXRKCB.findsortidx(x,  5, 7, 9) == 0
    @test IMEXRKCB.findsortidx(x,  9, 5, 6) == 0
end

@testset "Test functionality                     " begin
    # data
    N = 2
    ts = 0.0:0.1:1

    # create storage with one sample
    sol = IMEXRKCB.SolutionStorage(zeros(N))

    # push a cubic function of time
    for t in ts
        push!(sol, [t*t*t*j for j = 1:N], t)
    end

    # check interpolation
    out = zeros(N)
    for ti in 0.0:0.01:1.0
        out = sol(out, ti)
        @test out ≈ [ti*ti*ti, 2*ti*ti*ti]
    end

    # check out of bound range. Give a 0 idxcurr here, does not matter
    @test_throws ErrorException sol(out, -1.0)
    @test_throws ErrorException sol(out,  1.1)

    # storage must be writable
    @test_throws ErrorException push!(setrmode!(sol), zeros(N), 0)
    @test iswmode(setrmode!(sol)) == false
    @test iswmode(setwmode!(sol)) == true
    @test isrmode(setrmode!(sol)) == true
    @test isrmode(setwmode!(sol)) == false

    # test resetting
    reset!(sol)
    @test length(sol.ts) == 0
    @test length(sol.xs) == 0
    @test sol.idx == 0
end