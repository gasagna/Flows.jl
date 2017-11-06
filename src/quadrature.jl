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