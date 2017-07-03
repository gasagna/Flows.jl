export Monitor, reset!

"""
    m = Monitor((f1, f2, ..., fn), xq, [sizehint=100])

Construct an object that monitors and stores the value of functions 
`f1, f2, ... fn` along with time integration. The functions depend 
on the state `x` and optionally on the quadrature 
variable `q`. The argument `xq` is used so that the types of the outputs
`f1(xq), f2(xq), ..., fn(xq)` can be calculated, to allocate appropriate
storage. 

`Monitor` objects have two user-exposed fields:

    m.times
    m.samples

The first  contains the times instant at which samples have been stored. The
second is a tuple of vectors containing the samples associated to the functions, 
e.g., `m.samples[1]`, contains the samples of the first function.

If only one quantity is monitored, the associated function needs to be wrapped
in a one element tuple, e.g. as

    m = Monitor((f,), xq)

and the associated samples will be stored in `m.samples[1]`.

For functions that depend on the quadrature variable, the constructor is 
as follows

    m = Monitor((f,), (x, q))

i.e., state and quadrature variable instances are wrapped in a tuple, and 
the function `f` has signature

    function f(xq)
        # separate state and quadrature
        x, q = xq 
        # calculations
    end
"""
struct Monitor{F, S, N}
    fs::F
    samples::S
    times::Vector{Float64}
end

# Outer constructor
function Monitor(fs::NTuple{N, Base.Callable}, xq, sizehint::Int=100) where N
    types   = map(f->typeof(f(xq)), fs)
    samples = map(T->sizehint!(T[], sizehint), types)
    times   = sizehint!(Float64[], sizehint)
    Monitor{typeof(fs), typeof(samples), N}(fs, samples, times)
end

# push new samples to the monitor, at time t. This need to be a generated function
# to avoid allocations using a normal for loop. Try again with simpler code in
# future Julia versions...
@generated function Base.push!(m::Monitor{F, S, N}, t::Real, xq) where {F, S, N}
    expr = Expr(:block)
    for i = 1:N
        push!(expr.args, :(push!(m.samples[$i], m.fs[$i](xq))))
    end
    push!(expr.args, :(push!(m.times, t)))
    return expr
end

# Reset the data in a monitor.
function reset!(m::Monitor{F, S, N}, sizehint::Int=100) where {F, S, N}
    sizehint!(resize!(m.times, 0), sizehint)
    foreach(s->sizehint!(resize!(s, 0), sizehint), m.samples)
    return m
end