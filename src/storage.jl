export RAMStorage, times, samples

# ////// SOLUTION STORAGE //////
abstract type AbstractStorage{T, X} end

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

@inline Base.push!(rs::RAMStorage{T, X}, t::T, x::X) where {T, X} =
    (push!(rs.ts, t); push!(rs.xs, x); nothing)

times(  rs::RAMStorage) = rs.ts
samples(rs::RAMStorage) = rs.xs