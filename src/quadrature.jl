import Base: A_mul_B!

# This defines a type that wraps the main state vector `x::X` augmented
# with an appropriate object `q::Q` used for pure quadrature integration.
# Note that the type is immutable, hence if a scalar function need to be 
# integrated, its value is stored in a one-element vector
struct AugmentedState{X, Q}
    x::X
    q::Q
end

aug_state(x, q) = AugmentedState(x, q)
aug_state(x)    = x

# extract parts from augmented state object
_state(x::AugmentedState) = x.x
_quad(x::AugmentedState) = x.q
_state(x) = x
_quad(x) = x
_state_quad(x::AugmentedState) = (x.x, x.q)
_state_quad(x) = x

# Operations are broadcasted to both state and quadrature parts
@generated function Base.Broadcast.broadcast!(f, dest::AugmentedState, args...)
    quote
        $(Expr(:meta, :inline))
        broadcast!(f, _state(dest), map(_state, args)...)
        broadcast!(f,  _quad(dest), map(_quad,  args)...)
        dest
    end
end    

# Used in the definition of a scheme for proper allocation of storage. Note
# that `Y` cannot be a `Number`, because `similar` would raise an error. The
# state needs to be augmented with a mutable object!
Base.similar(z::AugmentedState) = AugmentedState(similar.(_state_quad(z))...)

# define this to broadcast the call to the subsystems, the main rhs of the
# ode, and the quadrature rule.
struct AugmentedSystem{G, A, Q}
    g::G # explicit term
    A::A # linear implicit term: can be Void if not provided
    q::Q # quadrature function: can be Void if not provided
end
aug_system(g, A, q) = AugmentedSystem(g, A, q)

function (aug::AugmentedSystem{G, A, Q})(t::Real, z::T, dzdt::T) where {G, A, Q, T} 
    aug.g(t, _state(z), _state(dzdt))
    Q <: Void || aug.q(t, _state(z), _quad(dzdt))
    return dzdt
end

function A_mul_B!(out::T, aug::AugmentedSystem{G, A}, z::T) where {G, A, T} 
    if A <: Void
        out .*= 0
    else
        A_mul_B!(_state(out), aug.A, _state(z))
    end
    return out
end

# Since we treat the quadrature fully explicitly, the solution of
# (I-cA)z = y for the quadrature part is simply z = y, because the
# component of A associated to this part is zero and the state 
# and quadrature parts are decoupled.
function ImcA!(aug::AugmentedSystem{G, A, Q}, c::Real, y::T, z::T) where {G, A, Q, T} 
    if A <: Void
        z .= y
    else
        ImcA!(aug.A, c, _state(y), _state(z))
        T<:AugmentedState && _quad(z) .= _quad(y)
    end
    return z
end