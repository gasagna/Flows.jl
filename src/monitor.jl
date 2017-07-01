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
struct Monitor{T, F} <: AbstractVector{Tuple{Float64, T}}
    f::F
    samples::Vector{T}
    times::Vector{Float64}
end

# Outer constructor
function Monitor(f, x, q...)
    # get the type of the output. This calls `f` once. Hope
    # it is not a problem....
    T = typeof(f(x, q...))
    Monitor{T, typeof(f)}(f, T[], Float64[])
end

Base.size(m::Monitor) = (length(m.times), )
Base.getindex(m::Monitor, i::Int) = (m.times[i], m.samples[i])

# push a new sample to the monitor, at time t
Base.push!(m::Monitor, t::Real, x, q...) = 
    (push!(m.samples, m.f(x, q...)); push!(m.times, t); nothing)

# if we have more monitors push to all of them
function Base.push!(ms::Tuple{Vararg{Monitor}}, t::Real, x, q...)
    for m in ms
        push!(m, t, x, q...)
    end
    nothing
end