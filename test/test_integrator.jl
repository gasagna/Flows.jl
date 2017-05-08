using Base.Test
using IMEXRKCB

# Define a custom type that satisfy the required interface. 
# Note that for subtypes of AbstractVector `A_mul_B!` and `ImcA!`
# already work out of the box.
struct foo{T} <: AbstractVector{T}
    data::Vector{T}
end
foo(data::AbstractVector{T}) where T = foo{T}(data)

Base.size(f::foo) = size(f.data)
Base.getindex( f::foo,      i::Int) =  f.data[i]
Base.setindex!(f::foo, val, i::Int) = (f.data[i] = val)
Base.similar(f::foo) = foo(similar(f.data))

@testset "forwmap!" begin
    # make system
    g(t, x, ẋ) = (ẋ .= -0.5.*x; ẋ)
    A = Diagonal([-0.5])

    # try on standard vector and on custom type
    for x in [[1.0], foo{Float64}([1.0])]

        # integration scheme
        for scheme in [IMEXRK3R2R(IMEXRKCB3e, x, false),
                       IMEXRK3R2R(IMEXRKCB3c, x, false),
                       IMEXRK4R3R(IMEXRKCB4,  x, false)]

            # forward map
            ϕ = forwmap!(g, A, 1.0, 0.01123, scheme)

            # check relative error
            @test ((ϕ(x) - exp(-1))/exp(-1))[1] < 4.95e-8
            @test ((ϕ(x) - exp(-2))/exp(-2))[1] < 4.95e-8
            @test ((ϕ(x) - exp(-3))/exp(-3))[1] < 4.95e-8
            @test ((ϕ(x) - exp(-4))/exp(-4))[1] < 4.95e-8
            @test ((ϕ(x) - exp(-5))/exp(-5))[1] < 4.95e-8

            # try giving a different input
            @test_throws MethodError ϕ([1])
        end
    end
end

@testset "time step" begin
    @test IMEXRKCB.next_Δt(0.0, 1.0, 0.1) == 0.1
    @test IMEXRKCB.next_Δt(0.0, 1.0, 1.1) == 1.0
end