# This defines a type that wraps the main state vector `x::X` augmented
# with an appropriate object `y::Y` used for pure quadrature integration.
# Note that the type is immutable, hence if a scalar function need to be 
# integrated, its value is stored in a one-element vector
struct AugmentedState{X, Y}
    x::X
    y::Y
end

aug_state(x, y) = AugmentedState(x, y)
aug_state(x)    = x

# extract state part from an augmented state
_state(x::AugmentedState) = x.x
_state(x) = x
_quad(x::AugmentedState) = x.y
_quad(x) = ()

# Used in the definition of a scheme for proper allocation of storage. Note
# that `Y` cannot be a `Number`, because `similar` would raise an error. The
# state needs to be augmented with a mutable object!
Base.similar(z::AugmentedState) = 
    AugmentedState(similar(_state(z)), similar(_quad(z)))

# define this to broadcast the call to the subsystems, the main rhs of the
# ode, and the quadrature rule.
struct AugmentedSystem{G, Q}
    g::G
    q::Q
end
aug_system(g, q) = AugmentedSystem(g, q)

(f::AugmentedSystem)(t::Real, z::AugmentedState, ż::AugmentedState) =
    (f.g(t, _state(z), _state(ż)); 
     f.q(t, _state(z), _quad(ż)))

# Broadcast calls to state part only
Base.A_mul_B!(out::AugmentedState, A, z::AugmentedState) = 
        (A_mul_B!(_state(out), A, _state(z)); _quad(out) .= 0)

# Since we treat the quadrature fully implicitly, the solution of
# (I-cA)z = y for the quadrature part is simply z = y, because the
# component of A associated to this part is zero and the state 
# and quadrature parts are decoupled.
ImcA!(A, c::Real, y::T, z::T) where T<:AugmentedState =
    (ImcA!(A, c, _state(y), _state(z)); _quad(z) .= _quad(y))