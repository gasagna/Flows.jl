@testset "linear system                          " begin

    # define linear system ẋ = x, but splitting the right hand side
    # into the explicit and implicit parts, to check both.
    g(t, x, ẋ) = (ẋ .= 0.5.*x; ẋ)
    A = Diagonal([0.5])

    #                                    impl       tableau     embed?  ord  bounds                embedded bounds
    for (impl, order, bnd, bnd_emb) in [(IMEXRK3R2R(IMEXRKCB2,  true),  2,  (0.025000, 0.027100), (0.017000, 0.020000)),
                                        (IMEXRK3R2R(IMEXRKCB3e, false), 3,  (0.006900, 0.008100), (0.000000, 0.000000)),
                                        (IMEXRK3R2R(IMEXRKCB3c, true),  3,  (0.007300, 0.008600), (0.023000, 0.027000)),
                                        (IMEXRK4R3R(IMEXRKCB4,  true),  4,  (0.000037, 0.000076), (0.001000, 0.001470))]

        # define scheme by implementation and storage type                                         
        scheme = IMEXRKScheme(impl, Float64[1.0])

        # ensure that the error decays with expected rate
        for Δt = [5^(-i) for i in linspace(1, 3, 5)]

            # step forward
            x = [1.0]
            IMEXRKCB.step!(scheme, IMEXRKCB.System(g, A, nothing), 0., Δt, x)

            # check error decays with expected power. The bounds bnd are used
            # to check whether the error decays at the expected rate, up
            # to a small variation, within the bounds.
            err = abs(x[1] - exp(Δt))
            # @show err/Δt^(order + 1)
            @test err/Δt^(order + 1) > bnd[1]
            @test err/Δt^(order + 1) < bnd[2]

            # test embedded scheme has order one less
            if isembedded(impl)
                x̂ = scheme.storage[end]
                err_emb = abs(x̂[1] - exp(Δt))
                # @show err_emb/Δt^order
                @test err_emb/Δt^order > bnd_emb[1]
                @test err_emb/Δt^order < bnd_emb[2]
            end
        end
    end
end