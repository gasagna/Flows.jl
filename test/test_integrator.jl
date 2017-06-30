# Define a custom type that satisfy the required interface. 
# Note that for subtypes of AbstractVector `A_mul_B!` and `ImcA!`
# already work out of the box.
struct foo{T}
    data::Vector{T}
end
foo(data::AbstractVector{T}) where T = foo{T}(data)

Base.eachindex(f::foo) = eachindex(f.data)
Base.similar(f::foo) = foo(similar(f.data))
Base.getindex( f::foo,      i::Int) =  f.data[i]
Base.setindex!(f::foo, val, i::Int) = (f.data[i] = val)

Base.A_mul_B!(out::foo{Float64}, A::Diagonal{Float64}, in::foo{Float64}) = 
    A_mul_B!(out.data, A, in.data)

@testset "integrator                             " begin
    # make system
    g(t, x, ẋ) = (for i in eachindex(x); ẋ[i] = -0.5*x[i]; end; ẋ)
    A = Diagonal([-0.5])

    # try on standard vector and on custom type
    for x in [[1.0], foo{Float64}([1.0])]

        # integration scheme
        for scheme in [IMEXRK3R2R(IMEXRKCB3e, false, x),
                       IMEXRK3R2R(IMEXRKCB3c, false, x),
                       IMEXRK4R3R(IMEXRKCB4,  false, x)]

            # forward map
            ϕ = integrator(g, A, scheme, 0.01123)

            # check relative error
            @test ((ϕ(x, 1.0)[1] - exp(-1))/exp(-1)) < 4.95e-8
            @test ((ϕ(x, 1.0)[1] - exp(-2))/exp(-2)) < 4.95e-8
            @test ((ϕ(x, 1.0)[1] - exp(-3))/exp(-3)) < 4.95e-8
            @test ((ϕ(x, 1.0)[1] - exp(-4))/exp(-4)) < 4.95e-8
            @test ((ϕ(x, 1.0)[1] - exp(-5))/exp(-5)) < 4.95e-8

            # try giving a different input
            @test_throws MethodError ϕ([1], 1.0)
        end
    end
end

@testset "time step                              " begin
    @test IMEXRKCB.next_Δt(0.0, 1.0, 0.1) == 0.1
    @test IMEXRKCB.next_Δt(0.0, 1.0, 1.1) == 1.0
end