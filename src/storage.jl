export setrmode!, setwmode!, SolutionStorage, iswmode, isrmode, reset!

# Storage for forward solution
abstract type AbstractStorage{X} end

# Store/Read solution along forward/backward problem
mutable struct SolutionStorage{X} <: AbstractStorage{X}
    xs::Vector{X}       # snapshots
    ts::Vector{Float64} # times
    tag::Symbol         # change mode of operation
    idx::Int            # index of last snapshot used
end

# Outer constructor makes the storage writeable by default
SolutionStorage(x::X) where {X} = SolutionStorage{X}(X[], Float64[], :w, 0)

# Add snapshots and time to the storage
function Base.push!(sol::SolutionStorage{X}, x::X, t::Real) where {X}
    sol.tag == :w || throw(ErrorException("storage must be in write mode. Use "))
    push!(sol.xs, x); push!(sol.ts, t); sol.idx += 1
end

# reset the storage
reset!(sol::SolutionStorage) = 
    (resize!(sol.ts, 0); resize!(sol.xs, 0); sol.idx =0; return sol)

# Forward and backward problems do different things to the storage
iswmode(sol::SolutionStorage) = sol.tag == :w
isrmode(sol::SolutionStorage) = sol.tag == :r
setwmode!(sol::SolutionStorage) = (sol.tag = :w; sol)
setrmode!(sol::SolutionStorage) = (sol.tag = :r; sol)

# Functor for 3rd order Lagrange interpolation
function (sol::SolutionStorage{X})(out::X, t::Real) where {X}
    # fix small round off error
    t < sol.ts[1]   && (t = t + 1e-10)
    t > sol.ts[end] && (t = t - 1e-10)
    # check if t is inbounds
    (t < sol.ts[1] || t > sol.ts[end]) && throw(error("selected time $t is out" *
                                     " of range [$(sol.ts[1]), $(sol.ts[end])]"))

    # Use current idx to speed up search. If not found (==0) use full length.
    idx = findsortidx(sol.ts, t, max(1, sol.idx-5), min(sol.idx+5, length(sol.ts)))
    idx == 0 && (idx = findsortidx(sol.ts, t))

    # boundary conditions need shifting of the stencil
    Δ = idx == 1                  ?  1 :
        idx == length(sol.ts)     ? -2 :
        idx == length(sol.ts) - 1 ? -1 : 0
    idx += Δ

    # get interpolation weights
    w1, w2, w3, w4 = cubicLagrangeWeights(sol.ts[idx-1],
                                          sol.ts[idx],
                                          sol.ts[idx+1],
                                          sol.ts[idx+2], t)
    # compute linear combination
    out .= w1.*sol.xs[idx-1] .+ w2.*sol.xs[idx] .+
           w3.*sol.xs[idx+1] .+ w4.*sol.xs[idx+2]

    # update idx to restart the next search close to where we left
    sol.idx = idx

    return out
end

@inline cubicLagrangeWeights(x0, x1, x2, x3, x) =
    ((x-x1)*(x-x2)*(x-x3)/((x0-x1)*(x0-x2)*(x0-x3)),
     (x-x0)*(x-x2)*(x-x3)/((x1-x0)*(x1-x2)*(x1-x3)),
     (x-x0)*(x-x1)*(x-x3)/((x2-x0)*(x2-x1)*(x2-x3)),
     (x-x0)*(x-x1)*(x-x2)/((x3-x0)*(x3-x1)*(x3-x2)))

# Returns the index i of sequence ts for which t > ts[i] and t < ts[i+1].
# This is an O(log n) function that expects a sorted sequence. Returns 0
# in case of something is going wrong, to signal that `t` does not fall
# between `ts[low]` and `ts[high]`.
@inline function findsortidx{T<:Integer}(ts::AbstractVector, x::Real, low::T, high::T)
    # error cases
    high < low && return 0
    x < ts[low] && return 0
    x > ts[high] && return 0
    # boundary cases
    x == ts[low] && return low
    x == ts[high] && return high
    # general exit cases
    low+one(low) == high && return low
    # split the interval
    mid = low + (high-low)>>1
    x > ts[mid] ? findsortidx(ts, x, mid, high) : findsortidx(ts, x, low, mid)
end
@inline findsortidx(ts, x) = findsortidx(ts, x, 1, length(ts))