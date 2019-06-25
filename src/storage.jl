export RAMStorage, times, samples, degree

# ////// SOLUTION STORAGE //////
"""
    AbstractStorage{T, X, DEG}

Abstract type for storage of solution data. The type parameters specify:
`T<:Real`    : the type used to store times
`X`          : the type used to store snapshots
`DEG`        : the degree of the interpolating lagrange polynomial
"""
abstract type AbstractStorage{T<:Real, X, DEG} end

degree(::AbstractStorage{T, X, DEG}) where {T, X, DEG} = DEG

# checkers
_isstorage(::Type{<:AbstractStorage}) = true
_isstorage(::Any) = false

# interface for subtypes
reset!(store::AbstractStorage, sizehint::Int) = error("not implemented")
Base.push!(store::AbstractStorage{T, X}, t::T, x::X) where {T, X} =
    error("not implemented")

times(  store::AbstractStorage) = error("not implemented")
samples(store::AbstractStorage) = error("not implemented")


# /// RAM storage ///
struct RAMStorage{T,
                  X,
                  DEG,
                  Vt<:AbstractVector{T},
                  Vx<:AbstractVector{X}} <: AbstractStorage{T, X, DEG}
        ts::Vt # Vector of times
        xs::Vx # Vector of snapshots
    period::T  # The period in case data is periodic. Defaults to `0.0`, so non periodic.
               # If the period is finite, it is assumed that the data is uniformly spaced
               # and that the last/first element is not repeated in the sequence

    # constructor from type
    function RAMStorage(::Type{X};
                   ttype::Type{T}=Float64,
                  degree::Int=3,
                  period::Real=0.0) where {T, X}
        degree ≥ 3 || throw(ArgumentError("polynomial degree must be ≥ 3, got $degree."))
        degree % 2 == 1 || throw(ArgumentError("polynomial degree must be odd, got $degree"))
        period ≥ 0 || throw(ArgumentError("period must be non-negative, got $period"))
        return new{T, X, degree, Vector{T}, Vector{X}}(T[], X[], period)
    end

    # constructor from object
    RAMStorage(::X; kwargs...) where {X} = RAMStorage(X; kwargs...)
end

@inline reset!(rs::RAMStorage, sizehint::Int=0) =
    (sizehint!(empty!(rs.ts), sizehint); sizehint!(empty!(rs.xs), sizehint); rs)

@inline Base.push!(rs::RAMStorage{T, X}, t::Real, x::X) where {T, X} =
    (push!(rs.ts, t); push!(rs.xs, x); nothing)

times(  rs::RAMStorage) = rs.ts
samples(rs::RAMStorage) = rs.xs

# /// LAGRANGIAN INTERPOLATION ///

"""
    _lagr_weights(t::Real, ts::NTuple{N, Real}) where {N}

Compute a `N`-tuple for the interpolation weights to interpolate
data at position `t` using data points available at positions `ts`,
with a lagrange polynomial of degree `N-1`.
"""
@generated _lagr_weights(t::Real, ts::NTuple{N, Real}, ::Val{0}) where {N} =
    :(Base.Cartesian.@ntuple $N j->_prod(t, ts, Val(j))/_prod(ts[j], ts, Val(j)))

"""
    _prod(t::T, ts::NTuple{N, T}, ::Val{SKIP}) where {N, T, SKIP}

Compute the product `(t - ts[1])*(t - ts[2])...(t - ts[N])` excluding the
factor `(t - ts[SKIP])`.
"""
@generated _prod(t::T, ts::NTuple{N, T}, ::Val{SKIP}) where {N, T, SKIP} =
    :(return *($([:(t - ts[$k]) for k in 1:N if k != SKIP]...)))

"""
    _lagr_interp(out::X, t::Real, ts::NTuple{N, Real},
                 xs::AbstractVector{X}, rng::NTuple{N, Int}) where {X, N}

Interpolate data points `(ts, xs[rng])` at location `t`, using a
lagrange polynomial of degree `N-1`, and overwrite the first argument `out`.
"""
function _lagr_interp(out::X,
                       t::Real,
                      ts::NTuple{N, Real},
                      xs::AbstractVector{X},
                     rng::NTuple{N, Int},
                   order::Val{ORD}) where {X, N, ORD}

    # get weights
    ws = _lagr_weights(t, ts, order)

    # compute linear combination and return output
    out .= ws[1].*xs[rng[1]]
    for i in 2:N
        out .+= ws[i].*xs[rng[i]]
    end

    return out
end

# accept being out of bounds by this much!
const TTOL = 1e-10

"""
    _normalise(t::Real, T::Real, isperiodic::Bool)

If `isperiodic` is `true` normalise `t` to be between `0` and `T`, else just return `t`.
"""
_normalise(t::Real, T::Real, isperiodic::Bool) =
    isperiodic ? (t + (1 + abs(t ÷ T))*T) % T : t

"""
    _interp_range(t::Real, ts::AbstractVector{<:Real}, ::Val{N}, period::Real) where {N}

Return an `N`-tuple of integer indices for the elements of `ts` that
participate in the interpolation at some point `t`. If the underlying
data is periodic, provide the `period`, with `Inf` to signal non-periodicity.
"""
function _interp_range(t::Real,
                      ts::AbstractVector{<:Real},
                        ::Val{N},
                  period::Real) where {N}
    #  we must have enough data
    length(ts) ≥ N ||
        throw(ArgumentError("input array length must be greater than N"))

    # determine if data is periodic
    isperiodic = iszero(period) ? false : true

    # If the data is periodic, we reset time `t` to be between `0` and `period`
    t = _normalise(t, period, isperiodic)

    # If the data is not periodic, check if `t` is in bounds. Note that although
    # time should never goes out of bounds in normal integration, it might be
    # possible in the runge kutta steps that we do so by a very small amount,
    # but generally only for values larger than the end of `ts`.
    if isperiodic == false
        isbetween(t, ts[1], ts[end]+TTOL) ||
            throw(ArgumentError("selected time $t is out of range [$(ts[1]), $(ts[end])]"))
    end

    # search index
    idx = searchsortedlast(ts, t)

    # boundary conditions in the non-periodic case need shifting of the stencil
    if isperiodic == false
        Δ = idx == 0              ?  N>>1     : # this can only happen if the above check fails
            idx == 1              ?  N>>1 - 1 :
            idx == length(ts)     ? -N>>1     :
            idx == length(ts) - 1 ? -N>>1 + 1 : 0
        idx += Δ
    end

    return ntuple(N) do j
        mod(idx - N>>1 + j - 1 + length(ts), length(ts)) + 1
    end
end

"""
    (store::RAMStorage{T, X, DEG})(out::X, t::Real) where {T, X, DEG}

Interpolate the storage data at time `t` and overwrite the first argument
`out`, using lagrangian interpolation of order `DEG`.
"""
function (store::RAMStorage{T, X, DEG})(out::X,
                                          t::Real,
                                      order::Val{ORD}=Val(0)) where {T,
                                                                     X,
                                                                     DEG,
                                                                     ORD}
    # Aliases. These should be lazy objects
    ts, xs = times(store), samples(store)

    # Obtain the indices of the elements that participate in the interpolation
    idxs = _interp_range(t, ts, Val(DEG+1), store.period)

    # construct a tuple of increasing times
    _ts = ntuple(DEG+1) do j
        ts[idxs[j]] + store.period
    end

    # call interp function
    return _lagr_interp(out, t, _ts, xs, idxs, order)
end