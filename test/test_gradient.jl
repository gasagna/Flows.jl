import LinearAlgebra: Diagonal, norm, dot
using Statistics

# ---------------------------------------------------------------------------- #
# NONLINEAR EQUATIONS. We arbitrarily split the equations into explicit and
# implicit components, for testing purposes. To do so, we integrate the 
# diagonal part implicitly, which does also not depend on the state. 
# If you set flag to zero, you do not need the linear part.
struct Lorenz
    flag::Int
    force::Real
end

@inline function (eq::Lorenz)(t::Real, u::V, dudt::V) where {V <: AbstractVector}
    x, y, z = u
    dudt[1] =   10 * (y - x)      - eq.flag*( - 10*x ) + eq.force
    dudt[2] =   28 *  x - y - x*z - eq.flag*( - y )
    dudt[3] = -8/3 * z + x*y      - eq.flag*( - 8/3*z )
    return dudt
end

# ---------------------------------------------------------------------------- #
# TANGENT EQUATIONS

struct LorenzTan
    flag::Int
end

(eq::LorenzTan)(t::Real, u::V, dudt::V,
    v::V, dvdt::V) where {V<:AbstractVector} = eq(t, u, v, dvdt)

function (eq::LorenzTan)(t::Real, u::V,
                                  v::V, dvdt::V) where {V<:AbstractVector}
    # extract components
    x′, y′, z′ = v
    x,  y,  z  = u

    dvdt[1] =  10 * (y′ - x′)        - eq.flag*( - 10*x′ ) + 1
    dvdt[2] =  (28-z)*x′ - y′ - x*z′ - eq.flag*( - y′    )
    dvdt[3] = -8/3*z′ + x*y′ + x′*y  - eq.flag*( - 8/3*z′)

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

    dvdt[1] =  -10*x⁺ - (z - 28)*y⁺ +    y*z⁺ - eq.flag*( - 10*x⁺ )
    dvdt[2] =   10*x⁺ -          y⁺ +    x*z⁺ - eq.flag*( - y⁺    )
    dvdt[3] =         -        x*y⁺ -  8/3*z⁺ - eq.flag*( - 8/3*z⁺) + 1
    
    return dvdt
end

# The diagonal term is integrated implicitly
const A = Diagonal([-10, -1, -8/3])

function trapz(xs, ys)
    N = length(xs)
    length(ys) == length(xs) || throw(ArgumentError("argh!!!"))
    I = (ys[2]+ys[1])*(xs[2]-xs[1])
    for k = 2:N-1
        I += (ys[k+1]+ys[k])*(xs[k+1]-xs[k])
    end
    return I/2
end

# @testset "test2                    " begin
#     T = 1
    
#     for dt = [1e-2, 1e-3, 1e-4]
#         x0 = Float64[9.1419853, 1.648665, 35.21793]
#         y0 = zeros(3)
#         w0 = zeros(3)

#         METHOD = RK4
#         flag = 0 

#         # finite difference
#         method = METHOD(x0, :NORMAL)
#         alpha = 1e-7
#         ϕp = flow(Lorenz(flag, +alpha), method, TimeStepConstant(dt))
#         ϕm = flow(Lorenz(flag, -alpha), method, TimeStepConstant(dt))

#         mon = Monitor(x0, x->x[3])
#         ϕp(copy(x0), (0, T), reset!(mon))
#         vp = trapz(times(mon), samples(mon))
#         ϕm(copy(x0), (0, T), reset!(mon))
#         vm = trapz(times(mon), samples(mon))
#         Jp_FD = (vp-vm)/2/alpha

#         # tangent
#         method = METHOD(couple(x0, y0), :NORMAL)
#         ψ = flow(couple(Lorenz(flag, 0), LorenzTan(flag)), method, TimeStepConstant(dt))
#         mon = Monitor(couple(x0, y0), x->x[2][3])
#         ψ(couple(copy(x0), y0), (0, T), reset!(mon))
#         Jp_TAN = trapz(times(mon), samples(mon))
#         # println(times(mon))

#         # adjoint
#         method = METHOD(x0, :NORMAL)
#         ϕ = flow(Lorenz(flag, 0), method, TimeStepConstant(dt))
#         cache = RAMStageCache(4, similar(x0))
#         ϕ(copy(x0), (0, T), cache)
        
#         method = METHOD(w0, :ADJ)
#         ψ_A = flow(LorenzAdj(flag), method, TimeStepFromCache())
#         mon = Monitor(w0, x->x[1])
#         ψ_A(w0, cache, reset!(mon))
#         Jp_ADJ = -trapz(times(mon), samples(mon))
#         # println(times(mon))
#         # println(samples(mon))

#         println(Jp_FD)
#         println(Jp_TAN)
#         println(Jp_ADJ) 
#         println( abs(Jp_ADJ - Jp_TAN)/abs(Jp_TAN)/dt^2 )
#         println()
#     end
# end

@testset "test                    " begin
    T = 1
        
    for dt = [1e-2, 1e-3, 1e-4]
        x0 = Float64[9.1419853, 1.648665, 35.21793]
        y0 = zeros(3)
        w0 = zeros(3)

        METHOD = CNRK2
        flag = 0
        IMPL = flag*A

        # finite difference
        method = METHOD(x0, :NORMAL)
        alpha = 1e-7
        ϕp = flow(Lorenz(flag, +alpha), IMPL, method, TimeStepConstant(dt))
        ϕm = flow(Lorenz(flag, -alpha), IMPL, method, TimeStepConstant(dt))

        mon = Monitor(x0, x->x[3])
        ϕp(copy(x0), (0, T), reset!(mon))
        vp = trapz(times(mon), samples(mon))
        ϕm(copy(x0), (0, T), reset!(mon))
        vm = trapz(times(mon), samples(mon))
        Jp_FD = (vp-vm)/2/alpha

        # tangent
        method = METHOD(couple(x0, y0), :NORMAL)
        ψ = flow(couple(Lorenz(flag, 0), LorenzTan(flag)), couple(IMPL, IMPL), method, TimeStepConstant(dt))
        mon = Monitor(couple(x0, y0), x->x[2][3])
        ψ(couple(copy(x0), y0), (0, T), reset!(mon))
        Jp_TAN = trapz(times(mon), samples(mon))
        # println(times(mon))

        # adjoint
        method = METHOD(x0, :NORMAL)
        ϕ = flow(Lorenz(flag, 0), IMPL, method, TimeStepConstant(dt))
        cache = RAMStageCache(2, similar(x0))
        ϕ(copy(x0), (0, T), cache)
        
        method = METHOD(w0, :ADJ)
        ψ_A = flow(LorenzAdj(flag), IMPL, method, TimeStepFromCache())
        mon = Monitor(w0, x->x[1])
        ψ_A(w0, cache, reset!(mon))
        Jp_ADJ = -trapz(times(mon), samples(mon))
        # println(times(mon))
        # println(samples(mon))

        println(Jp_FD)
        println(Jp_TAN)
        println(Jp_ADJ) 
        println( abs(Jp_ADJ - Jp_TAN)/abs(Jp_TAN)/dt^2 )
        println()
    end
end