using Base.Test
using Flows

# ---------------------------------------------------------------------------- #
# NONLINEAR EQUATIONS. We arbitrarily split the equations into explicit and
# implicit components, for testing purposes. To do so, we integrate the 
# diagonal part implicitly, which does also not depend on the state. 
# If you set flag to zero, you do not need the linear part.
struct Lorenz
    flag::Int
end

@inline function (eq::Lorenz)(t::Real, u::V, dudt::V) where {V <: AbstractVector}
    x, y, z = u
    @inbounds dudt[1] =   10 * (y - x)      - eq.flag*( - 10*x )
    @inbounds dudt[2] =   28 *  x - y - x*z - eq.flag*( - y )
    @inbounds dudt[3] = -8/3 * z + x*y      - eq.flag*( - 8/3*z )
    return dudt
end

# ---------------------------------------------------------------------------- #
# TANGENT EQUATIONS

struct LorenzTan
    flag::Int
end

# define two methods, with and without dudt
(eq::LorenzTan)(t::Real, u::V, dudt::V,
                         v::V, dvdt::V) where {V<:AbstractVector} =
    eq(t, u, v, dvdt)

function (eq::LorenzTan)(t::Real, u::V,
                                  v::V, dvdt::V) where {V<:AbstractVector}
    # extract components
    x′, y′, z′ = v
    x,  y,  z  = u

    @inbounds dvdt[1] =  10 * (y′ - x′)        - eq.flag*( - 10*x′ )
    @inbounds dvdt[2] =  (28-z)*x′ - y′ - x*z′ - eq.flag*( - y′    )
    @inbounds dvdt[3] = -8/3*z′ + x*y′ + x′*y  - eq.flag*( - 8/3*z′)

    return dvdt
end


# ---------------------------------------------------------------------------- #
# ADJOINT EQUATIONS

struct LorenzAdj
    flag::Int
end

function (eq::LorenzAdj)(t::Real, u::V, v::V, dvdt::V) where {V<:AbstractVector}
    # extract components
    x⁺, y⁺, z⁺ = v
    x,  y,  z  = u

    @inbounds dvdt[1] =  -10*x⁺ - (z - 28)*y⁺ +    y*z⁺ - eq.flag*( - 10*x⁺ )
    @inbounds dvdt[2] =   10*x⁺ -          y⁺ +    x*z⁺ - eq.flag*( - y⁺    )
    @inbounds dvdt[3] =         -        x*y⁺ -  8/3*z⁺ - eq.flag*( - 8/3*z⁺)
    
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
    nl    = RK4(x0, :NORMAL)
    l_t = RK4(x0, :TAN)
    l_adj = RK4(x0, :ADJ)

    # stage cache
    scache = RAMStageCache(4, x0)

    # system (without diagonal)
    sys_nl    = Flows.System(Lorenz(0),    nothing, nothing)
    sys_l_tan = Flows.System(LorenzTan(0), nothing, nothing)
    sys_l_adj = Flows.System(LorenzAdj(0), nothing, nothing)

    # execute step
    N = 5
    for i = 1:N
        Flows.step!(nl, sys_nl,    0, 1e-2, x0, scache)
    end

    y0 = Float64[1, 2, 3]
    for i = 1:N
        Flows.step!(l_t, sys_l_tan, 0, 1e-2, y0, scache.xs[i])
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
    # tb = @belapsed Flows.step!($l_t, $sys_l_tan, 0, 1e-2, $q1, $(scache.xs[1]))
    # println(ta/tb)
end

@testset "CB3R2R                                 " begin 
    # initial conditions
    x0 = Float64[15, 16, 20]

    for (nl, l_t, l_adj, NS) in [(CB3R2R2(x0,  :NORMAL),  CB3R2R2(x0, :TAN),  CB3R2R2(x0, :ADJ),  3),
                                 (CB3R2R3e(x0, :NORMAL), CB3R2R3e(x0, :TAN), CB3R2R3e(x0, :ADJ), 4),
                                 (CB3R2R3c(x0, :NORMAL), CB3R2R3c(x0, :TAN), CB3R2R3c(x0, :ADJ), 4)]

        # stage cache
        scache = RAMStageCache(NS, x0)

        # system
        sys_nl    = Flows.System(Lorenz(1),    A, nothing)
        sys_l_tan = Flows.System(LorenzTan(1), A, nothing)
        sys_l_adj = Flows.System(LorenzAdj(1), A, nothing)

        # execute step
        N = 50
        for i = 1:N
            Flows.step!(nl, sys_nl,    0, 1e-2, x0, scache)
        end

        y0 = Float64[1, 2, 3]
        for i = 1:N
            Flows.step!(l_t, sys_l_tan, 0, 1e-2, y0, scache.xs[i])
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
        # tb = @belapsed Flows.step!($l_t, $sys_l_tan, 0, 1e-2, $q1, $(scache.xs[1]))
        # println(ta/tb)
    end
end

# ---------------------------------------------------------------------------- #
# TEST LINEARISED STEP IS REALLY THE LINEARISATION OF THE NONLINEAR STEP
@testset "Complex step derivative                " begin

    # initial conditions using complex numbers
    x0 = zeros(Complex128, 3)

    # complex step
    ϵ = 1e-100

    for (nl, l_t, NS, _g_nl, _g_t, _A) in [(RK4(x0,      :NORMAL),      RK4(real.(x0), :TAN), 4, Lorenz(0), LorenzTan(0), nothing),
                                           (CB3R2R2(x0,  :NORMAL),  CB3R2R2(real.(x0), :TAN), 3, Lorenz(1), LorenzTan(1), A),
                                           (CB3R2R3e(x0, :NORMAL), CB3R2R3e(real.(x0), :TAN), 4, Lorenz(1), LorenzTan(1), A),
                                           (CB3R2R3c(x0, :NORMAL), CB3R2R3c(real.(x0), :TAN), 4, Lorenz(1), LorenzTan(1), A)]

        for i = 1:3
            x0 = [9.1419853, 1.648665, 35.21793] + im*[0.0, 0.0, 0.0]

            # stage cache
            scache = RAMStageCache(NS, x0)

            # system (without diagonal)
            sys_nl    = Flows.System(_g_nl, _A, nothing)
            sys_l_tan = Flows.System(_g_t,  _A, nothing)

            # go to attractor
            for j = 1:100
                Flows.step!(nl, sys_nl, 0, 1e-2, x0, nothing)
            end

            # reset perturbation
            x0 .= real.(x0)
            x0[i] += ϵ*im

            # number of steps
            N = 1000

            for j = 1:N
                Flows.step!(nl, sys_nl, 0, 1e-2, x0, scache)
            end

            # perturb z component only
            y0 = Float64[0, 0, 0]
            y0[i] += 1
            for j = 1:N
                Flows.step!(l_t, sys_l_tan, 0, 1e-2, y0, real.(scache.xs[j]))
            end

            @test norm(imag.(x0)./ϵ - y0)/norm(y0) < 5e-14
        end
    end
end

# ---------------------------------------------------------------------------- #
# TEST LINEARISED EQUATIONS API
@testset "linear api                             " begin

    # initial conditions
    x0 = Float64[9.1419853, 1.648665, 35.21793]

    for (nl, l_t, l_a, NS, _g_nl, _g_t, _g_a, _A) in [(RK4(x0,      :NORMAL),      RK4(x0, :TAN), RK4(x0,      :ADJ), 4, Lorenz(0), LorenzTan(0), LorenzAdj(0), nothing),
                                                      (CB3R2R2(x0,  :NORMAL),  CB3R2R2(x0, :TAN), CB3R2R2(x0,  :ADJ), 3, Lorenz(1), LorenzTan(1), LorenzAdj(1), A),
                                                      (CB3R2R3e(x0, :NORMAL), CB3R2R3e(x0, :TAN), CB3R2R3e(x0, :ADJ), 4, Lorenz(1), LorenzTan(1), LorenzAdj(1), A),
                                                      (CB3R2R3c(x0, :NORMAL), CB3R2R3c(x0, :TAN), CB3R2R3c(x0, :ADJ), 4, Lorenz(1), LorenzTan(1), LorenzAdj(1), A)]


        # stage cache
        scache = RAMStageCache(NS, x0)

        # non linear flow map
        ϕ = flow(_g_nl, _A, nl, TimeStepConstant(1e-2))

        # linearised propagator and adjoint
        ψ  = flow(_g_t, _A, l_t, TimeStepFromCache())
        ψ⁺ = flow(_g_a, _A, l_a, TimeStepFromCache())

        # propagate nonlinear operator
        ϕ(copy(x0), (0, 5), reset!(scache))

        # propagate linear operators forward/backward
        y0 = Float64[1, 2, 3]
        ψ(y0, scache)
        qT = Float64[4, 5, 6]
        ψ⁺(qT, scache)

        # test 
        a = dot(y0, [4, 5, 6])
        b = dot(qT, [1, 2, 3])
        @test abs(a - b)/a < 1e-14
    end
end