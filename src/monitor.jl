import DataStructures: searchsortedlast

export Monitor, reset!

# ///  Abstract type for all solution monitors ///
abstract type AbstractMonitor{T, X} end

# ///// UTILS //////
# whether t is between low and high
isbetween(t::Real, low::Real, high::Real) = (t ≥ low && t ≤ high)

_ismonitor(::Type{<:AbstractMonitor}) = true
_ismonitor(::Any) = false


# /// Monitor to save all time steps ///
mutable struct Monitor{T, X, S<:AbstractStorage{T, X}, F} <: AbstractMonitor{T, X}
          store::S                       # (time, samples) tuples
              f::F                       # action on what is begin pushed
       oneevery::Int                     # save every ... time steps
    savebetween::Tuple{Float64, Float64} # save only between these two times 
          count::Int                     # how many items we have in the store
    Monitor(store::S,
                f::F,
                oneevery::Int, 
                savebetween::Tuple{Real, Real}) where {T, X, S<:AbstractStorage{T, X}, F} =
        new{T, X, S, F}(store, f, oneevery, savebetween, 0)
end

# Provide a sample of what will be pushed
Monitor(x,
        f::Base.Callable=identity,
        store::S=RAMStorage(f(x));
        oneevery::Int=1,
        savebetween::Tuple{Real, Real}=(-Inf, Inf),
        sizehint::Int=0) where {S<:AbstractStorage} =
    Monitor(reset!(store, sizehint), f, oneevery, savebetween)

# Add sample and time to the storage
@inline function Base.push!(mon::Monitor, t::Real, x, force::Bool=false)
    if force == true || (mon.count % mon.oneevery == 0)
        if isbetween(t, mon.savebetween...)
            push!(mon.store, t, mon.f(x))
        end
    end
    mon.count += 1
    return nothing
end

# Reset storage
reset!(mon::Monitor, sizehint::Int=0) =
    (reset!(mon.store, sizehint); mon.count = 0; mon)

# get times or samples
times(mon::Monitor)   = times(mon.store)
samples(mon::Monitor) = samples(mon.store)