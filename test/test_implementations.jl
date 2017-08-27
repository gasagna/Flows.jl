@testset "isembedded                             " begin
    for embed in [true, false]
        for impl in [IMEXRK3R2R(IMEXRKCB2, embed), 
                     IMEXRK3R2R(IMEXRKCB3e, embed), 
                     IMEXRK3R2R(IMEXRKCB3c, embed),
                     IMEXRK4R3R(IMEXRKCB4,  embed)]
            @test isembedded(impl)         == embed
            @test isembedded(typeof(impl)) == embed
        end
    end
end

@testset "nstorage                               " begin
    @test nstorage(IMEXRK3R2R(IMEXRKCB2,  false)) == 3
    @test nstorage(IMEXRK3R2R(IMEXRKCB3e, false)) == 3
    @test nstorage(IMEXRK3R2R(IMEXRKCB3c, false)) == 3
    @test nstorage(IMEXRK4R3R(IMEXRKCB4,  false)) == 4
    @test nstorage(IMEXRK3R2R(IMEXRKCB2,  true))  == 4
    @test nstorage(IMEXRK3R2R(IMEXRKCB3e, true))  == 4
    @test nstorage(IMEXRK3R2R(IMEXRKCB3c, true))  == 4
    @test nstorage(IMEXRK4R3R(IMEXRKCB4,  true))  == 5
    @test nstorage(typeof(IMEXRK3R2R(IMEXRKCB2,  false))) == 3
    @test nstorage(typeof(IMEXRK3R2R(IMEXRKCB3e, false))) == 3
    @test nstorage(typeof(IMEXRK3R2R(IMEXRKCB3c, false))) == 3
    @test nstorage(typeof(IMEXRK4R3R(IMEXRKCB4,  false))) == 4
    @test nstorage(typeof(IMEXRK3R2R(IMEXRKCB2,  true)))  == 4
    @test nstorage(typeof(IMEXRK3R2R(IMEXRKCB3e, true)))  == 4
    @test nstorage(typeof(IMEXRK3R2R(IMEXRKCB3c, true)))  == 4
    @test nstorage(typeof(IMEXRK4R3R(IMEXRKCB4,  true)))  == 5
end