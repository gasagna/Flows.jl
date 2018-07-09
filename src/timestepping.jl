export AbstractTimeStepping,
       AbstractTimeStepFromHook,
       TimeStepConstant

abstract type AbstractTimeStepping end
abstract type AbstractTimeStepFromHook <: AbstractTimeStepping end

struct TimeStepConstant <: AbstractTimeStepping
    Δt::Float64
    function TimeStepConstant(Δt::Real)
        Δt > 0 || throw(ArgumentError("time step must be positive"))
        new(Float64(Δt))
    end
end