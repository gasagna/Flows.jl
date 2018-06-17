@testset "nstages                                " begin
    t = Flows.Tableau([1 2;
                          3 4],
                         [5, 6],
                         [7, 8],
                         [9, 0.])
    @test Flows.nstages(t) == 2
end

@testset "promote                                " begin
    t = Flows.Tableau([1 2;
                          3 4],
                         [5, 6],
                         [7, 8],
                         [9, 0.])
    @test typeof(t[:a, 1, 1]) == Float64
    @test typeof(t[:b, 1]   ) == Float64
    @test typeof(t[:e, 1]   ) == Float64
    @test typeof(t[:c, 1]   ) == Float64
end

@testset "convert                                " begin
    t = Flows.Tableau([1 2;
                          3 4],
                         [5, 6],
                         [7, 8],
                         [9, 0.])
    for T in [Float64, Float32, Rational{Int}]
        tf = convert(Flows.Tableau{T}, t)
        @test typeof(tf[:a, 1, 1]) == T
        @test typeof(tf[:b, 1]   ) == T
        @test typeof(tf[:e, 1]   ) == T
        @test typeof(tf[:c, 1]   ) == T
    end
end

@testset "indexing                               " begin
    t = Flows.Tableau([1 2;
                          3 4],
                         [5, 6],
                         [7, 8],
                         [9, 0])

    @test t[:a, 1, 1] == 1
    @test t[:a, 1, 2] == 2
    @test t[:a, 2, 1] == 3
    @test t[:a, 2, 2] == 4
    @test t[:b, 1   ] == 5
    @test t[:b, 2   ] == 6
    @test t[:e, 1   ] == 7
    @test t[:e, 2   ] == 8
    @test t[:c, 1   ] == 9
    @test t[:c, 2   ] == 0
end

@testset "IMEXTableau convert                    " begin
    for T in [Int128, Int64, Float32, Float64]
        c = one(T)
        ti = Flows.Tableau([1  2;
                               3  4],
                              [5, 6],
                              [7, 8],
                              [9, 0])

        te = Flows.Tableau([11  12;   #a
                               13  14],
                              [15, 16],  #b
                              [17, 18],  #e
                              [19,  c])   #c

        # constructor
        t  = Flows.IMEXTableau(ti, te)
        tf = convert(Flows.IMEXTableau{Float64}, t) 
        @test typeof(tf) == Flows.IMEXTableau{Float64}
    end
end

@testset "IMEXTableau promote                    " begin
    for T in [Int128, Int64, Float32, Float64]
        c = one(T)
        ti = Flows.Tableau([1  2;
                               3  4],
                              [5, 6],
                              [7, 8],
                              [9, 0])

        te = Flows.Tableau([11  12;   #a
                               13  14],
                              [15, 16],  #b
                              [17, 18],  #e
                              [19,  c])   #c
        # constructor
        t = Flows.IMEXTableau(ti, te)
        @test typeof(t) == Flows.IMEXTableau{T}
    end
end

@testset "IMEXTableau indexing                   " begin
    ti = Flows.Tableau([1  2;
                           3  4],
                          [5, 6],
                          [7, 8],
                          [9, 0])

    te = Flows.Tableau([11  12;   #a
                           13  14],
                          [15, 16],  #b
                          [17, 18],  #e
                          [19, 10])  #c

    # constructor
    t = Flows.IMEXTableau(ti, te)

    # implicit
    @test t[:aᴵ, 1, 1] == 1
    @test t[:aᴵ, 1, 2] == 2
    @test t[:aᴵ, 2, 1] == 3
    @test t[:aᴵ, 2, 2] == 4
    @test t[:bᴵ, 1]    == 5
    @test t[:bᴵ, 2]    == 6
    @test t[:eᴵ, 1]    == 7
    @test t[:eᴵ, 2]    == 8
    @test t[:cᴵ, 1]    == 9
    @test t[:cᴵ, 2]    == 0

    # explicit
    @test t[:aᴱ, 1, 1] == 11
    @test t[:aᴱ, 1, 2] == 12
    @test t[:aᴱ, 2, 1] == 13
    @test t[:aᴱ, 2, 2] == 14
    @test t[:bᴱ, 1]    == 15
    @test t[:bᴱ, 2]    == 16
    @test t[:eᴱ, 1]    == 17
    @test t[:eᴱ, 2]    == 18
    @test t[:cᴱ, 1]    == 19
    @test t[:cᴱ, 2]    == 10
end