using Base.Test
using IMEXRK

@testset "linear system" begin


    g(t, x, ẋ) = (ẋ .= 0.5.*x; ẋ)
    A = Diagonal([0.5])

    # make scheme and storage of appropriate type
    for scheme in [IMEXRK3R2R(IMEXRKCB3e, [1.0]),
                   IMEXRK3R2R(IMEXRKCB3c, [1.0])]

        for Δt = [5^(-i) for i in linspace(1, 4.5, 10)]

            # step forward
            x = [1.0]
            IMEXRK.step!(scheme, g, A, 0., Δt, x)

            # check error decays with 4th power
            err = abs(x[1] - exp(Δt))
            @test err/Δt^4 > 0.0068
            @test err/Δt^4 < 0.00852
        end
    end
end