@testset "lossy                                  " begin
    rng = Flows.Steps(0, 1, 0.4)

    expected = [(0.0, 0.4), (0.4, 0.4), (0.8, 0.2)]
    for (i, (t, dt)) in enumerate(rng)
        @test  t ≈ expected[i][1]
        @test dt ≈ expected[i][2]
    end

    # test interface used in the integration routines
    for rep = 1:100
        T  = 100 + rand()
        dt =       rand()
        rng = Flows.Steps(0, T, dt)
        for (i, (t, dt)) in enumerate(rng)
            @test t + dt ≤ T
        end
    end
end

@testset "lossless                               " begin
    rng = Flows.Steps(0, 0.8, 0.4)

    expected = [(0.0, 0.4), (0.4, 0.4)]
    for (i, (t, dt)) in enumerate(rng)
        @test  t ≈ expected[i][1]
        @test dt ≈ expected[i][2]
    end

    # test interface used in the integration routines
    T  = 10.0
    dt = 1e-3
    rng = Flows.Steps(0, T, dt)
    for (i, (t, dt)) in enumerate(rng)
        @test t + dt ≤ T
    end
end

@testset "backwards checks                       " begin
    @test_throws ArgumentError Flows.Steps(1, 0, -0.1)
end

@testset "backwards lossy                        " begin
    rng = Flows.Steps(1, 0, 0.4)
    expected = [(1.0, -0.4), (0.6, -0.4), (0.2, -0.2)]
    for (i, (t, dt)) in enumerate(rng)
        @test  t ≈ expected[i][1]
        @test dt ≈ expected[i][2]
    end

    # test interface used in the integration routines
    for rep = 1:100
        T  = 100 + rand()
        dt =       rand()
        rng = Flows.Steps(T, 0, dt)
        for (i, (t, dt)) in enumerate(rng)
            @test t + dt ≥ 0
        end
    end

end

@testset "backwards lossless                     " begin
    rng = Flows.Steps(0.8, 0, 0.4)
    
    expected = [(0.8, -0.4), (0.4, -0.4)]
    for (i, (t, dt)) in enumerate(rng)
        @test  t ≈ expected[i][1]
        @test dt ≈ expected[i][2]
    end
    
    # test interface used in the integration routines
    T  = 10.0
    dt = 1e-3
    rng = Flows.Steps(T, 0, dt)
    for (i, (t, dt)) in enumerate(rng)
        @test t + dt  ≥ 0
    end
end
