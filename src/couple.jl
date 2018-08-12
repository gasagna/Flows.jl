# This is basically a 2-tuple. We could use tuples, obtaining most of
# the functionality except `similar` and `copy`, which we would need to
# overload. That would be type piracy, so we create our implementation.
struct Couple{A, B}
    a::A
    b::B
end

# constructors
coupled(a, b) = Couple(a, b)

# extract parts
@inline Base.first(ab::Couple) = ab.a
@inline Base.last( ab::Couple) = ab.b

# Operations are broadcasted to both parts
@generated function Base.Broadcast.broadcast!(f, dest::Couple, args...)
    quote
        $(Expr(:meta, :inline))
        broadcast!(f, first(dest), map(first, args)...)
        broadcast!(f,  last(dest), map(last,  args)...)
        return dest
    end
end

Base.similar(ab::Couple) = coupled(similar(first(ab)), similar(last(ab)))
Base.copy(   ab::Couple) = coupled(   copy(first(ab)),    copy(last(ab)))