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
struct AugmentedSystem{G, Q}
    g::G
    q::Q
end
aug_system(g, q) = AugmentedSystem(g, q)

# treat the quadrature equation fully explicitly
(f::AugmentedSystem)(t::Real, z::AugmentedState, ż::AugmentedState) =
    (f.g(t, _state(z), _state(ż)); 
     f.q(t, _state(z), _quad(ż)))

# Broadcast calls to state part only
Base.A_mul_B!(out::T, A, z::T) where {T<:AugmentedState} =
    (A_mul_B!(_state(out), A, _state(z)); _quad(out) .= 0)

# fix ambiguity warning due to method for Diagonal in imca.jl
Base.A_mul_B!(out::T, A::Diagonal, z::T) where {T<:AugmentedState} =
    (A_mul_B!(_state(out), A, _state(z)); _quad(out) .= 0)

# Since we treat the quadrature fully explicitly, the solution of
# (I-cA)z = y for the quadrature part is simply z = y, because the
# component of A associated to this part is zero and the state 
# and quadrature parts are decoupled.
ImcA!(A, c::Real, y::T, z::T) where T<:AugmentedState =
    (ImcA!(A, c, _state(y), _state(z)); _quad(z) .= _quad(y))

# fix ambiguity warning due to method for Diagonal in imca.jl
ImcA!(A::Diagonal, c::Real, y::T, z::T) where T<:AugmentedState =
    (ImcA!(A, c, _state(y), _state(z)); _quad(z) .= _quad(y))    