export Monitor

"""
    m = Monitor(f, x[, q])

Construct an object that monitors and stores the value of a function 
`f` depending on the state `x` and optionally on the quadrature 
variable `q`. A monitor object has two user-exposed fields:

    m.times
    m.samples

which contains the times instant at which samples have been stored 
and the samples themselves, respectively.
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