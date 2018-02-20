# Define a custom type that satisfy the required interface.
# Note that for subtypes of AbstractVector `A_mul_B!` and `ImcA!`
# already work out of the box.
struct foo{T}
    data::Vector{T}
end
foo(data::AbstractVector{T}) where T = foo{T}(data)

Base.size(f::foo) = size(f.data)
Base.similar(f::foo) = foo(similar(f.data))
Base.copy(f::foo) = foo(copy(f.data))
Base.eltype(::foo{T}) where {T} = T
Base.getindex( f::foo,      i::Int) =  f.data[i]
Base.setindex!(f::foo, val, i::Int) = (f.data[i] = val)
Base.unsafe_get(f::foo) = f.data

@generated function Base.Broadcast.broadcast!(f, dest::foo{T}, args::Vararg{Any, N}) where {T, N}
    quote
        broadcast!(f, unsafe_get(dest), map(unsafe_get, args)...)
        # note we return dest, and not the output of the above
        # broadcast!, which would be unsafe_get(dest), i.e. the
        # private field data in this case.
        return dest
    end
end

@testset "integrator                             " begin
    # make system
    g(t, x, ẋ) = (ẋ .= .-0.5.*x; ẋ)
    A = Diagonal([-0.5])

    # try on standard vector and on custom type
    for x in [Float64[1.0], foo{Float64}([1.0])]

    # integration scheme
        for (method, err) in [(IMEXMethod(:CB2_3R2R,  copy(x)), 2e-5),
                              (IMEXMethod(:CB3e_3R2R, copy(x)), 5e-8),
                              (IMEXMethod(:CB3c_3R2R, copy(x)), 6e-8),
                              (IMEXMethod(:CB4_4R3R,  copy(x)), 3e-12)]

            # forward map
            ϕ = integrator(g, A, method, 0.01123)

            # check relative error, for a few repetitions of the integration
            @test abs(ϕ(copy(x), (0, 1))[1] - exp(-1))/exp(-1) < err
            @test abs(ϕ(copy(x), (0, 2))[1] - exp(-2))/exp(-2) < err
            @test abs(ϕ(copy(x), (0, 3))[1] - exp(-3))/exp(-3) < err
            @test abs(ϕ(copy(x), (0, 4))[1] - exp(-4))/exp(-4) < err
            @test abs(ϕ(copy(x), (0, 5))[1] - exp(-5))/exp(-5) < err

            # try giving a different input
            @test_throws MethodError ϕ([1], (0, 1))
        end
    end
end