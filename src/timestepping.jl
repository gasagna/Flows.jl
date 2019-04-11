export AbstractTimeStepping,
       AbstractTimeStepFromHook,
       TimeStepConstant,
       TimeStepFromCache,
       TimeStepFromStorage

# ---------------------------------------------------------------------------- #
# The mother of all time stepping schemes
abstract type AbstractTimeStepping end

# ---------------------------------------------------------------------------- #
# When the time step depends on the solution
abstract type AbstractTimeStepFromHook <: AbstractTimeStepping end

# ---------------------------------------------------------------------------- #
# Constant time stepping
mutable struct TimeStepConstant <: AbstractTimeStepping
    Δt::Float64
    function TimeStepConstant(Δt::Real)
        Δt > 0 || throw(ArgumentError("time step must be positive"))
        new(Float64(Δt))
    end
end

# ---------------------------------------------------------------------------- #
# Provide time stepping based on a stage cache from a nonlinear solution
struct TimeStepFromCache <: AbstractTimeStepping end

# ---------------------------------------------------------------------------- #
# Provide time stepping based on a storage from a nonlinear solution
struct TimeStepFromStorage <: AbstractTimeStepping 
    Δt::Float64
    function TimeStepFromStorage(Δt::Real)
        Δt > 0 || throw(ArgumentError("time step must be positive"))
        new(Float64(Δt))
    end
end