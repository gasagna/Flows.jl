export AbstractTimeStepping,
       AbstractTimeStepFromHook,
       TimeStepConstant,
       TimeStepFromCache

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