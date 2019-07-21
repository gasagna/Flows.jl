export Coupled, couple, couplecopy

"""
    Coupled{N, ARGS<:NTuple{N, Any}} <: AbstractVector{Any}

Couple together `N` objects and store them internally in a Julia `ntuple` of size `N`.

`Coupled` objects have two purposes in this package. 

  * First, it is sometimes necessary to integrate several dynamical systems in a 
coupled fashion, such as, for instance, the integration of variational equations
along with a nonlinear problem. Another example is the integration of a dynamical 
system with the addition of a quadrature equation. `Flow` objects [Flow](@ref) can be
fed with `Coupled` objects and know how to coordinate the time stepping and manage 
dependencies. 

  * Second, in these cases one also needs to construct an integration method object 
that preallocates a state vector coupling the states of each dynamical system.

Coupled` objects are immutable, but their elements can be mutable, allowing the 
content to be modified. The individual elements of a `Coupled` object can be 
accessed using the standard indexing notation [`Base.getindex`](@ref), but the elements
cannot be changed, i.e. there is no `Base.setindex!` defined for `Coupled` objects.
"""
struct Coupled{N, ARGS<:NTuple{N, Any}} <: AbstractVector{Any}
    args::ARGS
    Coupled(args::ARGS) where {N, ARGS<:NTuple{N, Any}} = new{N, ARGS}(args)
end

"""
    couple(args...)

Create a `Coupled` object from the sequence of arguments `args`. 

# Examples
```jldoctest
# one dimensional linear system
f(t, x, dxdt) = (dxdt .= x; dxdt)

# quadrature equation to compute the integral of the state x
q(t, x, dxdt, q, dqdt) = (dydt[1] .= sin(x[1]); dqdt)

# define the coupled flow operator
F = flow(couple(f, g),
        RK4(couple(zeros(1), zeros(1))), 
        TimeStepConstant(1e-3))

# define an initial condition
xq = couple([1.0], [0.0])

# march forward. xq[2] will contain an approximation of the
# integral of sin(exp(x)) from 0 to 1, (see wolframalpha!)
@assert abs(F(xq, (0, 1)) - 0.87495) < 1e-5
```
"""
couple(args...) = Coupled(args)

"""
    couplecopy(N::Int, x)
    
Couple `N` copies of `x` together, created using `Base.deepcopy`.

# Examples
```jldoctest
x = couplecopy(3, [1, 2, 3])
@assert x[2] == [1, 2, 3]
```
"""
couplecopy(N::Int, x) = couple(ntuple(i->deepcopy(x), N)...)

"""
    Base.getindex(x::Coupled, i::Int)

Return the `i`-th element of a `Coupled` object. 

# Examples
```jldoctest
x = couple([23], [34])
@assert x[2] == [34]
```
"""
Base.@propagate_inbounds function Base.getindex(x::Coupled{N}, i::Int) where {N}
    @boundscheck 1 ≤ i ≤ N || throw(BoundsError())
    @inbounds val = x.args[i]
    return val
end

"""
    Base.similar(x::Coupled)

Call `Base.similar` on the elements of `x` and couple them together in a new object.
"""
Base.similar(x::Coupled{N}) where {N} = couple(ntuple(i->similar(x[i]), N)...)

"""
    Base.copy(x::Coupled)

Create a copy of `x` by coupling copies of its elements in a new object.
"""
Base.copy(x::Coupled{N})    where {N} = couple(ntuple(i->copy(x[i]), N)...)

"""
    Base.size(::Coupled{N}) where {N}

Return `(N, )`, the size of a `Coupled{N}` object, coupling `N` elements together.
"""
@inline Base.size(::Coupled{N}) where {N} = (N, )


# ~ BROADCASTING ~

# DEV NOTES: this code makes coupled objects support the dot broadcasting notation
# by forwarding the broadcasting to each of the objects held together. Broadcasting
# of algebraic functions is meant to work only for operands that are either numbers
# or coupled objects. 
const CoupledStyle = Broadcast.ArrayStyle{Coupled}
Base.BroadcastStyle(::Type{<:Coupled}) = CoupledStyle()

# The approach used here is taken from Chris Rackaukas' MultiScaleArrays package. It
# works directly with Julia's `Broadcast.Broadcasted` objects and extracts recursively
# the components of a `Coupled` object, applying the fused function to each of them.
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

# This creates a new `Broadcasted` object with the fused function and the `i`-th 
# elements of each `Coupled` object involved in the expression, as well any numbers.
@inline unpack(bc::Broadcast.Broadcasted, v::Val) =
    Broadcast.Broadcasted(bc.f, _unpack(bc.args, v))

# Here we need Val for type stability
@inline unpack(x::Coupled, ::Val{i}) where {i} = getfield(x.args, i)
@inline unpack(x::Number, i) = x

# The first item is `unpack` and not `_unpack` when we have nested `Broadcasted` objects.
@inline _unpack(args::Tuple, i) = (unpack(args[1], i), _unpack(Base.tail(args), i)...)
@inline _unpack(args::Tuple{Any}, i) = (unpack(args[1], i),)
@inline _unpack(args::Tuple{}, i) = ()