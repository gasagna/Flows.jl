import IMEXRKCB: Tableau, IMEXTableau, nstages

@testset "convert_tuple" begin
    typeof(IMEXRKCB.convert_tuple(Float64, 1)) == Float64
    typeof(IMEXRKCB.convert_tuple(Float64, (1, 1))) == Float64
    typeof(IMEXRKCB.convert_tuple(Float64, ((1, 1), (1, 1)))) == Float64
end

# conversion/promotion
@testset "To Float64" begin
    t = Tableau(((1, 2),
                 (3, 4)),
                 (5, 6),
                 (7, 8),
                 (9, 0.0))
    @test typeof(typeof(t)[Val{:a}, 1, 1]) == Float64
    @test typeof(typeof(t)[Val{:b}, 1])    == Float64
    @test typeof(typeof(t)[Val{:b̂}, 1])    == Float64
    @test typeof(typeof(t)[Val{:c}, 1])    == Float64
end

# standard Tableau indexing
let
    t = typeof(Tableau(((1, 2),   #a
                        (3, 4)),
                        (5, 6),   #b
                        (7, 8),   #b̂
                        (9, 0)))  #c

    @test nstages(t) == 2

    @test t[Val{:a}, 1, 1] == 1
    @test t[Val{:a}, 1, 2] == 2
    @test t[Val{:a}, 2, 1] == 3
    @test t[Val{:a}, 2, 2] == 4
    @test t[Val{:b}, 1   ] == 5
    @test t[Val{:b}, 2   ] == 6
    @test t[Val{:b̂}, 1   ] == 7
    @test t[Val{:b̂}, 2   ] == 8
    @test t[Val{:c}, 1   ] == 9
    @test t[Val{:c}, 2   ] == 0
end


# IMEX tableau indexing
let
    tᴵ = Tableau(((1, 2),  #a
                  (3, 4)),
                  (5, 6),  #b
                  (7, 8),  #b̂
                  (9, 0))  #c

    tᴱ = Tableau(((11, 12),  #a
                  (13, 14)),
                  (15, 16),  #b
                  (17, 18),  #b̂
                  (19, 10))  #c

    # constructor
    t = typeof(IMEXTableau(tᴵ, tᴱ))

    # implicit
    @test t[Val{:aᴵ}, 1, 1] == 1
    @test t[Val{:aᴵ}, 1, 2] == 2
    @test t[Val{:aᴵ}, 2, 1] == 3
    @test t[Val{:aᴵ}, 2, 2] == 4
    @test t[Val{:bᴵ}, 1   ] == 5
    @test t[Val{:bᴵ}, 2   ] == 6
    @test t[Val{:b̂ᴵ}, 1   ] == 7
    @test t[Val{:b̂ᴵ}, 2   ] == 8
    @test t[Val{:cᴵ}, 1   ] == 9
    @test t[Val{:cᴵ}, 2   ] == 0

    # explicit
    @test t[Val{:aᴱ}, 1, 1] == 11
    @test t[Val{:aᴱ}, 1, 2] == 12
    @test t[Val{:aᴱ}, 2, 1] == 13
    @test t[Val{:aᴱ}, 2, 2] == 14
    @test t[Val{:bᴱ}, 1   ] == 15
    @test t[Val{:bᴱ}, 2   ] == 16
    @test t[Val{:b̂ᴱ}, 1   ] == 17
    @test t[Val{:b̂ᴱ}, 2   ] == 18
    @test t[Val{:cᴱ}, 1   ] == 19
    @test t[Val{:cᴱ}, 2   ] == 10
end