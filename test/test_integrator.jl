# Define a custom type that satisfy the required interface. 
# Note that for subtypes of AbstractVector `A_mul_B!` and `ImcA!`
# already work out of the box.
struct foo{T}
    data::Vector{T}
end
foo(data::AbstractVector{T}) where T = foo{T}(data)

Base.size(f::foo) = size(f.data)
Base.similar(f::foo) = foo(similar(f.data))
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

        # integration scheme
    for (impl, err) in   [(IMEXRK3R2R(IMEXRKCB2,  false), 2e-5),
                          (IMEXRK3R2R(IMEXRKCB3e, false), 5e-8),
                          (IMEXRK3R2R(IMEXRKCB3c, false), 6e-8),
                          (IMEXRK4R3R(IMEXRKCB4,  false), 3e-12)]
    
        # try on standard vector and on custom type
        for x in [[1.0], foo{Float64}([1.0])]

            # construct scheme
            scheme = IMEXRKScheme(impl, x)                       

            # forward map
            ϕ = integrator(g, A, scheme, 0.01123)

            # check relative error, for a few repetitions of the integration
            @test abs(ϕ(x, (0, 1))[1] - exp(-1))/exp(-1) < err
            @test abs(ϕ(x, (0, 1))[1] - exp(-2))/exp(-2) < err
            @test abs(ϕ(x, (0, 1))[1] - exp(-3))/exp(-3) < err
            @test abs(ϕ(x, (0, 1))[1] - exp(-4))/exp(-4) < err
            @test abs(ϕ(x, (0, 1))[1] - exp(-5))/exp(-5) < err

            # try giving a different input
            @test_throws MethodError ϕ([1], (0, 1))
        end
    end
end