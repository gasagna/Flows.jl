export RAMStageCache

# ---------------------------------------------------------------------------- #
# Abstract type for all stage caches. The first parameters is the number
# of cached intermediate stage values, and the second is their type, e.g.
# some subtype of an AbstractArray.
abstract type AbstractStageCache{N, X} end

# ---------------------------------------------------------------------------- #
# Do not cache the stages
struct NoCache{N, X} <: AbstractStageCache{N, X} end
_iscache(::Any) = false
_iscache(::AbstractCache) = true

# ---------------------------------------------------------------------------- #
# Stage cache where all stages are cached in RAM.
struct RAMStageCache{N, X} <: AbstractStageCache{N, X}
     ts::Vector{Float64}
    Δts::Vector{Float64}
     xs::Vector{NTuple{N, X}}
end

@inline Base.push!(ss::RAMStageCache{N, X},
                    t::Real,
                   Δt::Real,
                    x::NTuple{N, X}) where {N, X} =
    (push!(ss.ts, t); push!(ss.Δts, Δt); push!(ss.xs, x); nothing)
