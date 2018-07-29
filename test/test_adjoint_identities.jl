using Base.Test
using Flows

# NonLinear equations
@inline function Lorenz(t::Real, u::V, dudt::V) where {V <: AbstractVector}
    x, y, z = u
    @inbounds dudt[1] =   10 * (y - x)
    @inbounds dudt[2] =   28 *  x - y - x*z
    @inbounds dudt[3] = -8/3 * z + x*y
    return dudt
end

# tangent equations
function LorenzTan(t::Real, u::V, v::V, dvdt::V) where {V<:AbstractVector}
    # extract components
    x′, y′, z′ = v
    x,  y,  z  = u

    @inbounds dvdt[1] =  10 * (y′ - x′)
    @inbounds dvdt[2] =  (28-z)*x′ - y′ - x*z′
    @inbounds dvdt[3] = -8/3*z′ + x*y′ + x′*y

    return dvdt
end

# adjoint equations
function LorenzAdj(t::Real, u::V, v::V, dvdt::V) where {V<:AbstractVector}
    # extract components
    x⁺, y⁺, z⁺ = v
    x,  y,  z  = u

    @inbounds dvdt[1] =  -10*x⁺ - (z - 28)*y⁺ +    y*z⁺
    @inbounds dvdt[2] =   10*x⁺ -          y⁺ +    x*z⁺
    @inbounds dvdt[3] =         -        x*y⁺ -  8/3*z⁺
    
    return dvdt
end

@testset "RK4                                    " begin 

    # initial conditions
    x0 = Float64[1, 1, 2]

    # stage cache
    scache = RAMStageCache(4, x0)

    # methods
    nl_yes = RK4(x0, true)
    nl_no  = RK4(x0, false)
    l_tan  = RK4_Tan(x0)
    l_adj  = RK4_Adj(x0)

    # system
    sys_nl    = Flows.System(Lorenz,    nothing, nothing)
    sys_l_tan = Flows.System(LorenzTan, nothing, nothing)
    sys_l_adj = Flows.System(LorenzAdj, nothing, nothing)

    # execute step
    N = 1000
    for i = 1:N
        Flows.step!(nl_yes, sys_nl,    0, 1e-2, x0, scache)
    end

    y0 = Float64[1, 2, 3]
    for i = 1:N
        Flows.step!(l_tan, sys_l_tan, 0, 1e-2, y0, scache.xs[i])
    end

    q1 = Float64[4, 5, 7]
    for i = N:-1:1
        Flows.step!(l_adj, sys_l_adj, 0, 1e-2, q1, scache.xs[i])
    end

    a = dot(y0, [4, 5, 7])
    b = dot(q1, [1, 2, 3])

    @test abs(a-b)/abs(a) < 1e-14
end