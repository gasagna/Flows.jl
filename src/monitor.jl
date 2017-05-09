export Monitor

struct Monitor{T, F} <: AbstractVector{Tuple{Float64, T}}
    f::F
    samples::Vector{T}
    times::Vector{Float64}
end

Base.size(m::Monitor) = (length(m.times), )
Base.getindex(m::Monitor, i::Int) = (m.times[i], m.samples[i])

function Monitor(x, f)
    # get the type of the output
    T = typeof(f(x))
    Monitor{T, typeof(f)}(f, T[], Float64[])
end

# push a new sample to the monitor, at time t
Base.push!(m::Monitor, x, t::Real) = 
    (push!(m.samples, m.f(x)); push!(m.times, t); nothing)

# if we have more monitors push to all of them
function Base.push!(ms::Tuple{Vararg{Monitor}}, x, t::Real)
    for m in ms
        push!(m, x, t)
    end
    nothing
end