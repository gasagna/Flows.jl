using Base.Test
using IMEXRKCB

@testset "n monitors                             " begin
    m = Monitor((x->x[1], x->x[1]^2), [1.0])
    
    push!(m, 0.0, [0.0])
    push!(m, 0.1, [1.0])
    push!(m, 0.2, [2.0])
    push!(m, 0.3, [3.0])

    @test m.times      == [0.0, 0.1, 0.2, 0.3]
    @test m.samples[1] == [0.0, 1.0, 2.0, 3.0]
    @test m.samples[2] == [0.0, 1.0, 4.0, 9.0]
end

@testset "allocation                             " begin
    g(t, x, ẋ) = (ẋ[1] = x[1]; ẋ)
    A = Diagonal([0.0])

    # integration scheme
    scheme = IMEXRK3R2R(IMEXRKCB3e, false, [0.0])

    # monitors
    m = Monitor((x->x[1], x->x[1]^2), [1.0])

    # forward map
    ϕ = integrator(g, A, scheme, 0.1)
    
    # initial condition
    x₀ = [0.0]

    # warm up
    ϕ(x₀, 1, m) 

    # does not allocate because we do not grow the arrays in the monitor
    @test (@allocated ϕ(x₀, 1, m)) == 0
end