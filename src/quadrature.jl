# This defines a type that wraps the main state vector `x::X` augmented
# with an appropriate object `q::Q` used for pure quadrature integration.
# Note that the type is immutable, hence if a scalar function need to be
# integrated, its value is stored in a one-element vector
struct Augmented{X, Q}
    x::X
    q::Q
end

augment(x, q) = Augmented(x, q)
augment(x)    = x

# extract parts from augmented state object
@inline _state(x::Augmented) = x.x
@inline  _quad(x::Augmented) = x.q
@inline _state(x) = x
@inline  _quad(x) = x

# Operations are broadcasted to both state and quadrature parts
@generated function Base.Broadcast.broadcast!(f, dest::Augmented, args...)
    quote
        $(Expr(:meta, :inline))
        broadcast!(f, _state(dest), map(_state, args)...)
        broadcast!(f,  _quad(dest), map(_quad,  args)...)
        return dest
    end
end

# Used in the definition of a scheme for proper allocation of storage. Note
# that `Y` cannot be a `Number`, because `similar` would raise an error. The
# state needs to be augmented with a mutable object!
Base.similar(z::Augmented) = Augmented(similar(_state(z)), similar(_quad(z)))