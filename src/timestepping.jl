export AbstractTimeStepping,
       AbstractTimeStepFromHook,
       TimeStepConstant,
       TimeStepFromCache,
       TimeStepFromStorage

# The mother of all time stepping schemes
abstract type AbstractTimeStepping end

"""
    AbstractTimeStepFromHook

Abstract type for time stepping schemes where the time step is determined at runtime.

See [`Flows.jl Adaptive time stepping`](@ref)
"""
abstract type AbstractTimeStepFromHook <: AbstractTimeStepping end

"""
    TimeStepConstant(Δt::Real)

Specify that integration should be performed with constant time step `Δt`.
"""
struct TimeStepConstant <: AbstractTimeStepping
    Δt::Float64
    function TimeStepConstant(Δt::Real)
        Δt > 0 || throw(ArgumentError("time step must be positive"))
        new(Float64(Δt))
    end
end

# Provide time stepping based on a stage cache from a nonlinear solution
struct TimeStepFromCache <: AbstractTimeStepping end

"""
    TimeStepFromStorage(Δt::Real)

Specify that integration should be performed with constant time step `Δt`, and 
that an [`AbstractStorage`](@ref) will be required.
"""
struct TimeStepFromStorage <: AbstractTimeStepping 
    Δt::Float64
    function TimeStepFromStorage(Δt::Real)
        Δt > 0 || throw(ArgumentError("time step must be positive"))
        new(Float64(Δt))
    end
end