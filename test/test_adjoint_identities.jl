using Base.Test
using Flows

# ---------------------------------------------------------------------------- #
# NONLINEAR EQUATIONS. We arbitrarily split the equations into explicit and
# implicit components, for testing purposes. To do so, we integrate the 
# diagonal part implicitly, which does also not depend on the state.

struct Lorenz end

@inline function (::Lorenz)(t::Real, u::V, dudt::V) where {V <: AbstractVector}
    x, y, z = u
    @inbounds dudt[1] =   10 * (y - x)      - ( - 10*x )
    @inbounds dudt[2] =   28 *  x - y - x*z - ( - y )
    @inbounds dudt[3] = -8/3 * z + x*y      - ( - 8/3*z )
    return dudt
end

# ---------------------------------------------------------------------------- #
# TANGENT EQUATIONS

struct LorenzTan end

function (::LorenzTan)(t::Real, u::V, v::V, dvdt::V) where {V<:AbstractVector}
    # extract components
    x′, y′, z′ = v
    x,  y,  z  = u

    @inbounds dvdt[1] =  10 * (y′ - x′)        - ( - 10*x′ )
    @inbounds dvdt[2] =  (28-z)*x′ - y′ - x*z′ - ( - y′    )
    @inbounds dvdt[3] = -8/3*z′ + x*y′ + x′*y  - ( - 8/3*z′)

    return dvdt
end


# ---------------------------------------------------------------------------- #
# ADJOINT EQUATIONS

struct LorenzAdj end

function (::LorenzAdj)(t::Real, u::V, v::V, dvdt::V) where {V<:AbstractVector}
    # extract components
    x⁺, y⁺, z⁺ = v
    x,  y,  z  = u

    @inbounds dvdt[1] =  -10*x⁺ - (z - 28)*y⁺ +    y*z⁺ - ( - 10*x⁺ )
    @inbounds dvdt[2] =   10*x⁺ -          y⁺ +    x*z⁺ - ( - y⁺    )
    @inbounds dvdt[3] =         -        x*y⁺ -  8/3*z⁺ - ( - 8/3*z⁺)
    
    return dvdt
end

# ---------------------------------------------------------------------------- #
# The diagonal term is integrated implicitly
const A = Diagonal([-10, -1, -8/3])

# ---------------------------------------------------------------------------- #
# TEST DISCRETE ADJOINT
@testset "RK4                                    " begin 

    # initial conditions
    x0 = Float64[1, 1, 2]

    # methods
    nl    = RK4(x0, :NL)
    l_tan = RK4(x0, :TAN)
    l_adj = RK4(x0, :ADJ)

    # stage cache
    scache = RAMStageCache(4, x0)

    # system (without diagonal)
    sys_nl    = Flows.System(Lorenz(),    nothing, nothing)
    sys_l_tan = Flows.System(LorenzTan(), nothing, nothing)
    sys_l_adj = Flows.System(LorenzAdj(), nothing, nothing)

    # execute step
    N = 5
    for i = 1:N
        Flows.step!(nl, sys_nl,    0, 1e-2, x0, scache)
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

    # these take the same time
    # ta = @belapsed Flows.step!($l_adj, $sys_l_adj, 0, 1e-2, $q1, $(scache.xs[1]))
    # tb = @belapsed Flows.step!($l_tan, $sys_l_tan, 0, 1e-2, $q1, $(scache.xs[1]))
    # println(ta/tb)
end

@testset "CB3R2R                                 " begin 
    # initial conditions
    x0 = Float64[1, 1, 2]

    # test 3e only
    nl    = CB3R2R3e(x0, :NL)
    l_tan = CB3R2R3e(x0, :TAN)
    l_adj = CB3R2R3e(x0, :ADJ)

    # stage cache
    scache = RAMStageCache(4, x0)

    # system
    sys_nl    = Flows.System(Lorenz(),    A, nothing)
    sys_l_tan = Flows.System(LorenzTan(), A, nothing)
    sys_l_adj = Flows.System(LorenzAdj(), A, nothing)

    # execute step
    N = 50
    for i = 1:N
        Flows.step!(nl, sys_nl,    0, 1e-2, x0, scache)
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

    # the adjoint code is ~30% slower here
    # ta = @belapsed Flows.step!($l_adj, $sys_l_adj, 0, 1e-2, $q1, $(scache.xs[1]))
    # tb = @belapsed Flows.step!($l_tan, $sys_l_tan, 0, 1e-2, $q1, $(scache.xs[1]))
    # println(ta/tb)
end