export RAMStorage, times, samples

# ////// SOLUTION STORAGE //////
abstract type AbstractStorage{T, X} end

# interface for subtypes
Base.sizehint!(store::AbstractStorage, hint::Int) = error("not implemented")
Base.resize!(store::AbstractStorage, size::Int) = error("not implemented")
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

@inline Base.sizehint!(rs::RAMStorage, hint::Int) = 
    (sizehint!(rs.ts, hint); sizehint!(rs.xs, hint); rs)
@inline Base.resize!(rs::RAMStorage, size::Int) = 
	(resize!(rs.ts, size); resize!(rs.xs, size); rs)
@inline Base.push!(rs::RAMStorage{T, X}, t::T, x::X) where {T, X} = 
    (push!(rs.ts, t); push!(rs.xs, x); nothing)

times(rs::RAMStorage) = rs.ts
samples(rs::RAMStorage) = rs.xs

# /// Lazy disk storage ///