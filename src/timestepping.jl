export AbstractTimeStepping, AbstractTimeStepFromHook

abstract type AbstractTimeStepping end
abstract type AbstractTimeStepFromHook <: AbstractTimeStepping end 

mutable struct TimeStepConstant <: AbstractTimeStepping
    Δt::Float64
    function ConstantTimeStep(Δt::Real)
        Δt > 0 || throw(ArgumentError("time step must be positive"))
        new(Float64(Δt))
    end
end


# change dt: replace with getproperty
set_Δt!(stepping::TimeStepConstant, Δt::Real) =
    (Δt < 0 && throw(ArgumentError("time step must be positive"));
     stepping.Δt = Δt; nothing)
