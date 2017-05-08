@testset "linear system" begin

    # define linear system
    g(t, x, ẋ) = (ẋ .= 0.5.*x; ẋ)
    A = Diagonal([0.5])

    #                                      impl       tableau    storage embed?  ord  bounds                embedded bounds
    for (scheme, order, bnd, bnd_emb) in [(IMEXRK3R2R(IMEXRKCB3e, [1.0], false), 3,  (0.006900, 0.008100), (0.000000, 0.000000)),
                                          (IMEXRK3R2R(IMEXRKCB3c, [1.0], true),  3,  (0.007300, 0.008600), (0.023000, 0.027000)),
                                          (IMEXRK4R3R(IMEXRKCB4,  [1.0], true),  4,  (0.000037, 0.000076), (0.001000, 0.001470))]

        # ensure that the error decays with expected rate
        for Δt = [5^(-i) for i in linspace(1, 3, 5)]

            # step forward
            x = [1.0]
            IMEXRKCB.step!(scheme, g, A, 0., Δt, x)

            # check error decays with expected power
            err = abs(x[1] - exp(Δt))
            @test err/Δt^(order + 1) > bnd[1]
            @test err/Δt^(order + 1) < bnd[2]

            # test embedded scheme has order one less
            if isembedded(scheme)
                x̂ = scheme.storage[end]
                err_emb = abs(x̂[1] - exp(Δt))
                @test err_emb/Δt^order > bnd_emb[1]
                @test err_emb/Δt^order < bnd_emb[2]
            end
        end
    end
end