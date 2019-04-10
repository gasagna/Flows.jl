export RAMStorage, times, samples

# ////// SOLUTION STORAGE //////
abstract type AbstractStorage{T, X} end

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
                  Vt<:AbstractVector{T},
                  Vx<:AbstractVector{X}} <: AbstractStorage{T, X}
    ts::Vt
    xs::Vx
    RAMStorage{T, X}() where {T, X} = new{T, X, Vector{T}, Vector{X}}(T[], X[])
end

@inline reset!(rs::RAMStorage, sizehint::Int=0) = 
    (sizehint!(empty!(rs.ts), sizehint); sizehint!(empty!(rs.xs), sizehint); rs)

@inline Base.push!(rs::RAMStorage{T, X}, t::Real, x::X) where {T, X} =
    (push!(rs.ts, t); push!(rs.xs, x); nothing)

times(  rs::RAMStorage) = rs.ts
samples(rs::RAMStorage) = rs.xs


# /// THIRD ORDER LAGRANGIAN INTERPOLATION ///

# get interpolation weights for the function
@inline weights(t, t0, t1, t2, t3, ::Val{0}) =
    ((t-t1)*(t-t2)*(t-t3)/((t0-t1)*(t0-t2)*(t0-t3)),
     (t-t0)*(t-t2)*(t-t3)/((t1-t0)*(t1-t2)*(t1-t3)),
     (t-t0)*(t-t1)*(t-t3)/((t2-t0)*(t2-t1)*(t2-t3)),
     (t-t0)*(t-t1)*(t-t2)/((t3-t0)*(t3-t1)*(t3-t2)))

# and its time derivative
@inline weights(t, t0, t1, t2, t3, ::Val{1}) =
    (((t-t1)*(t-t2) + (t-t1)*(t-t3) + (t-t2)*(t-t3))/((t0-t1)*(t0-t2)*(t0-t3)),
     ((t-t0)*(t-t2) + (t-t0)*(t-t3) + (t-t2)*(t-t3))/((t1-t0)*(t1-t2)*(t1-t3)),
     ((t-t1)*(t-t0) + (t-t1)*(t-t3) + (t-t0)*(t-t3))/((t2-t0)*(t2-t1)*(t2-t3)),
     ((t-t1)*(t-t2) + (t-t1)*(t-t0) + (t-t2)*(t-t0))/((t3-t0)*(t3-t1)*(t3-t2)))

#  The parameter deg indicates whether we want to interpolate 
# the function, deg=Val(0), or its time derivative, deg=Val(1) 
# (I can't remember why we need that).
function _lagr_3_interp(out::X,
                          t::Real,
                         x0::X,    x1::X,    x2::X,    x3::X,
                         t0::Real, t1::Real, t2::Real, t3::Real, deg::Val) where {X}
    # checks
    isbetween(t, min(t0, t3) - 1e-10, max(t0, t3)+1e-10) ||
        error("selected time is out of range")

    # get weights
    w0, w1, w2, w3 = weights(t, t0, t1, t2, t3, deg)

    # compute linear combination and return output
    out .= w0.*x0 .+ w1.*x1 .+ w2.*x2 .+ w3.*x3

    return out
end

function (mon::RAMStorage{T, X})(out::X, t::Real, deg::Val=Val(0)) where {T, X}
    # Aliases. These should be lazy objects
    ts, xs = times(mon), samples(mon)

    # check if t is inbounds. Note that although time never goes out the bounds of
    # the span, it might be possible in the runge kutta steps that
    # we do so by a very small amount, so we take care of that here.
    isbetween(t,  min(ts[1], ts[end]) - 1e-10, max(ts[1], ts[end]) + 1e-10) ||
        error("selected time is out of range")

    # search current index
    idx = searchsortedlast(ts, t)

    # boundary conditions need shifting of the stencil
    Δ = idx == 0              ?  2 : # this can only happen if the above check fails
        idx == 1              ?  1 :
        idx == length(ts)     ? -2 :
        idx == length(ts) - 1 ? -1 : 0
    idx += Δ

    # call interp function
    return _lagr_3_interp(out, t, xs[idx-1], xs[idx], xs[idx+1], xs[idx+2],
                                  ts[idx-1], ts[idx], ts[idx+1], ts[idx+2], deg)
end