# ---------------------------------------------------------------------------- #
# NONLINEAR EQUATIONS. We arbitrarily split the equations into explicit and
# implicit components, for testing purposes. To do so, we integrate the 
# diagonal part implicitly, which does also not depend on the state. 
# If you set flag to zero, you do not need the linear part.

struct Lorenz
    flag::Int
    force::Float64
    Lorenz(flag::Int, force::Real=0.0) = new(flag, force)
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
    force::Float64
    LorenzTan(flag::Int, force::Real=0.0) = new(flag, force)
end

(eq::LorenzTan)(t::Real, u::V, dudt::V,
    v::V, dvdt::V) where {V<:AbstractVector} = eq(t, u, v, dvdt)

function (eq::LorenzTan)(t::Real, u::V,
                                  v::V, dvdt::V) where {V<:AbstractVector}
    # extract components
    x′, y′, z′ = v
    x,  y,  z  = u

    dvdt[1] =  10 * (y′ - x′)        - eq.flag*( - 10*x′ ) + eq.force
    dvdt[2] =  (28-z)*x′ - y′ - x*z′ - eq.flag*( - y′    )
    dvdt[3] = -8/3*z′ + x*y′ + x′*y  - eq.flag*( - 8/3*z′)

    return dvdt
end


# ---------------------------------------------------------------------------- #
# ADJOINT EQUATIONS

struct LorenzAdj
    flag::Int
    force::Float64
    LorenzAdj(flag::Int, force::Real=0.0) = new(flag, force)
end

function (eq::LorenzAdj)(t::Real, u::V, v::V, dvdt::V) where {V<:AbstractVector}
    # extract components
    x⁺, y⁺, z⁺ = v
    x,  y,  z  = u

    dvdt[1] =  -10*x⁺ + (28 - z)*y⁺ +    y*z⁺ - eq.flag*( - 10*x⁺ )
    dvdt[2] =   10*x⁺ -          y⁺ +    x*z⁺ - eq.flag*( - y⁺    )
    dvdt[3] =         -        x*y⁺ -  8/3*z⁺ - eq.flag*( - 8/3*z⁺) 
    
    dvdt[3] += eq.force
    return 
end

# The diagonal term is integrated implicitly
const A = Diagonal([-10, -1, -8/3])