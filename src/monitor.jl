export Monitor, reset!

struct Monitor{X, V<:AbstractVector, F, O}
    ts::Vector{Float64} # times
    xs::V               # samples
     f::F               # act on what is begin pushed (copy, identity, other func)
    function Monitor{X, V, F, O}(ts::Vector, xs::V, f::F, sizehint::Int) where {X, V<:AbstractVector, F, O}
        length(ts) == length(xs) == 0 || error("input arrays must be zero")
        new{X, V, F, O}(sizehint!(ts, sizehint), sizehint!(xs, sizehint), f)
    end
end

# Provide a sample of what will be pushed
function Monitor(x::X, f::Base.Callable=identity, order::Int=3, sizehint::Int=100) where {X}
    T = typeof(f(x))
    Monitor{X, Vector{T}, typeof(f), order}(Float64[], T[], f, sizehint)
end

# Add snapshots and time to the storage
Base.push!(mon::Monitor{X}, t::Real, x::X) where {X} =
    (push!(mon.xs, mon.f(x)); push!(mon.ts, t))

# Reset storage
reset!(mon::Monitor, sizehint::Int=100) =
    (sizehint!(resize!(mon.ts, 0), sizehint);
     sizehint!(resize!(mon.xs, 0), sizehint); return mon)

# ~~~ Interpolation ~~~

# Third order interpolation
function (interp::Monitor{X, V, F, 3})(out::X, t::Real) where {X, V, F}
    # aliases
    xs, ts = interp.xs, interp.ts

    # fix small round off error
    t < ts[1]   && (t = t + 1e-10)
    t > ts[end] && (t = t - 1e-10)

    # check if t is inbounds
    t ≥ ts[1] && t ≤ ts[end] || throw(error("selected time $t is out" *
                                     " of range [$(ts[1]), $(ts[end])]"))

    # search current index
    idx = searchsortedlast(ts, t)

    # boundary conditions need shifting of the stencil
    Δ = idx == 1              ?  1 :
        idx == length(ts)     ? -2 :
        idx == length(ts) - 1 ? -1 : 0
    idx += Δ

    # get interpolation weights
    w1, w2, w3, w4 = weights(ts[idx-1],
                             ts[idx],
                             ts[idx+1],
                             ts[idx+2], t)

    # compute linear combination and return
    return out .= w1.*xs[idx-1] .+ w2.*xs[idx] .+ w3.*xs[idx+1] .+ w4.*xs[idx+2]
end

# Third order interpolation weights
@inline weights(t0, t1, t2, t3, t) =
    ((t-t1)*(t-t2)*(t-t3)/((t0-t1)*(t0-t2)*(t0-t3)),
     (t-t0)*(t-t2)*(t-t3)/((t1-t0)*(t1-t2)*(t1-t3)),
     (t-t0)*(t-t1)*(t-t3)/((t2-t0)*(t2-t1)*(t2-t3)),
     (t-t0)*(t-t1)*(t-t2)/((t3-t0)*(t3-t1)*(t3-t2)))