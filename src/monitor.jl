using DataStructures

export Monitor, reset!

# ///  Abstract type for all solution monitors ///
abstract type AbstractMonitor{T, X} end

# ///// UTILS //////
# whether t is between low and high
isbetween(t::Real, low::Real, high::Real) = (t ≥ low && t ≤ high)



# /// Monitor to save all time steps ///
mutable struct Monitor{T, X, O, S<:AbstractStorage{T, X}, F} <: AbstractMonitor{T, X}
       store::S  # (time, samples) tuples
           f::F  # action on what is begin pushed
    oneevery::Int
       count::Int
    Monitor{O}(store::S, 
               f::F, 
               oneevery::Int) where {T, X, O, S<:AbstractStorage{T, X}, F} =
        new{T, X, O, S, F}(store, f, oneevery, 0)
end

# Provide a sample of what will be pushed
Monitor(x,
        f::Base.Callable=identity,
        store::S = RAMStorage{Float64, typeof(f(x))}();
        oneevery::Int=1,
        order::Int=3,
        sizehint::Int=0) where {S<:AbstractStorage} =
    Monitor{order}(sizehint!(store, sizehint), f, oneevery)

# Add sample and time to the storage
@inline function Base.push!(mon::Monitor, t::Real, x)
    if mon.count % mon.oneevery == 0
        push!(mon.store, t, mon.f(x))
    end
    mon.count += 1
    return nothing
end

# Reset storage
reset!(mon::Monitor, sizehint::Int=0) =
    (sizehint!(resize!(mon.store, 0), sizehint); mon.count = 0; mon)

# get times or samples
times(mon::Monitor)   = times(mon.store)
samples(mon::Monitor) = samples(mon.store)



# /// Interpolation ///

# get interpolation weights
@inline weights(t, t0, t1, t2, t3, ::Val{0}) =
    ((t-t1)*(t-t2)*(t-t3)/((t0-t1)*(t0-t2)*(t0-t3)),
     (t-t0)*(t-t2)*(t-t3)/((t1-t0)*(t1-t2)*(t1-t3)),
     (t-t0)*(t-t1)*(t-t3)/((t2-t0)*(t2-t1)*(t2-t3)),
     (t-t0)*(t-t1)*(t-t2)/((t3-t0)*(t3-t1)*(t3-t2)))

@inline weights(t, t0, t1, t2, t3, ::Val{1}) =
    (((t-t1)*(t-t2) + (t-t1)*(t-t3) + (t-t2)*(t-t3))/((t0-t1)*(t0-t2)*(t0-t3)),
     ((t-t0)*(t-t2) + (t-t0)*(t-t3) + (t-t2)*(t-t3))/((t1-t0)*(t1-t2)*(t1-t3)),
     ((t-t1)*(t-t0) + (t-t1)*(t-t3) + (t-t0)*(t-t3))/((t2-t0)*(t2-t1)*(t2-t3)),
     ((t-t1)*(t-t2) + (t-t1)*(t-t0) + (t-t2)*(t-t0))/((t3-t0)*(t3-t1)*(t3-t2)))

# Third order Lagrangian interpolation
function lagrinterp(out::X,
                    t::Real,
                    x0::X,    x1::X,    x2::X,    x3::X,
                    t0::Real, t1::Real, t2::Real, t3::Real, deg::Val) where {X}
    # checks
    isbetween(t, min(t0, t3) - 1e-10, max(t0, t3)+1e-10) ||
        error("selected time is out of range")

    # get weights
    w0, w1, w2, w3 = weights(t, t0, t1, t2, t3, deg)

    # compute linear combination and return
    out .= w0.*x0 .+ w1.*x1 .+ w2.*x2 .+ w3.*x3

    return out
end

function (mon::Monitor{T, X, 3})(out::X, t::Real, deg::Val=Val{0}()) where {T, X}
    # Aliases. These should be lazy objects
    ts, xs = times(mon), samples(mon)

    # check if t is inbounds
    isbetween(t,  min(ts[1], ts[end]) - 1e-10, max(ts[1], ts[end])+1e-10) ||
        error("selected time is out of range")

    # search current index
    idx = searchsortedlast(ts, t)

    # boundary conditions need shifting of the stencil
    Δ = idx == 1              ?  1 :
        idx == length(ts)     ? -2 :
        idx == length(ts) - 1 ? -1 : 0
    idx += Δ

    # call interp function
    return lagrinterp(out, t, xs[idx-1], xs[idx], xs[idx+1], xs[idx+2],
                              ts[idx-1], ts[idx], ts[idx+1], ts[idx+2], deg)
end