export Coupled, couple, couplecopy

# This is basically an N-tuple. We could use Base.Tuple{Any, Any}, obtaining 
# most of the functionality except `similar` and `copy`, which we would need to
# overload. That would be type piracy, so we create our implementation. The 
# only hacky bit is how to define the broadcast behaviour for the Coupled object.
struct Coupled{N, ARGS<:NTuple{N, Any}} <: AbstractVector{Any}
    args::ARGS
    Coupled(args::ARGS) where {N, ARGS<:NTuple{N, Any}} = new{N, ARGS}(args)
end

# Constructor
couple(args...) = Coupled(args)

# couple `N` copies of `x`
couplecopy(N::Int, x) = couple(ntuple(i->deepcopy(x), N)...)

# Extract parts using indexing (this is read-only)
Base.@propagate_inbounds function Base.getindex(x::Coupled, i::Int)
    @boundscheck _checkbounds(x, i)
    @inbounds val = x.args[i]
    return val
end

_checkbounds(::Coupled{N}, i) where {N} = 
    (1 ≤ i ≤ N || throw(BoundsError()); nothing)

# Array interface
Base.similar(x::Coupled{N}) where {N} = couple(ntuple(i->similar(x[i]), N)...)
Base.copy(x::Coupled{N})    where {N} = couple(ntuple(i->copy(x[i]), N)...)
@inline Base.size(::Coupled{N})   where {N} = (N, )
@inline Base.length(::Coupled{N}) where {N} = N


# ~ BROADCASTING ~
const CoupledStyle = Broadcast.ArrayStyle{Coupled}
Base.BroadcastStyle(::Type{<:Coupled}) = CoupledStyle()

@generated function Base.copyto!(dest::Coupled{N},
                                   bc::Broadcast.Broadcasted{CoupledStyle}) where {N}
    quote
        $(Expr(:meta, :inline))
        Base.Cartesian.@nexprs $N i -> begin
            @inbounds copyto!(getfield(dest.args, i), unpack(bc, Val(i)))
        end
        return dest
    end
end

@inline unpack(bc::Broadcast.Broadcasted, v::Val) =
    Broadcast.Broadcasted(bc.f, _unpack(bc.args, v))

# we need Val for type stability
@inline unpack(x::Coupled, ::Val{i}) where {i} = getfield(x.args, i)
@inline unpack(x::Any, i) = x

@inline _unpack(args::Tuple, i) = (unpack(args[1], i), _unpack(Base.tail(args), i)...)
@inline _unpack(args::Tuple{Any}, i) = (unpack(args[1], i),)
@inline _unpack(args::Tuple{}, i) = ()