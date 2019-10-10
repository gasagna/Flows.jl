export ContinuousMode, DiscreteMode

# tag objects to specify the integration mode
abstract type AbstractMode end

struct NormalMode <: AbstractMode end

struct ContinuousMode{ISADJOINT} <: AbstractMode
    ContinuousMode(isadjoint::Bool=false) = new{isadjoint}()
end

struct DiscreteMode{ISADJOINT} <: AbstractMode
    DiscreteMode(isadjoint::Bool=false) = new{isadjoint}()
end

isadjoint(::Type{NormalMode}) = false
isadjoint(::Type{ContinuousMode{ISADJOINT}}) where {ISADJOINT} = ISADJOINT
isadjoint(::Type{DiscreteMode{ISADJOINT}})   where {ISADJOINT} = ISADJOINT

# Type parameters:
# - X           : the state type
# - MODE        : integration mode: normal, linear continuous, linear discrete
# - NS          : number of internal stages that are required to be saved to then 
#                 integrate linearised equations in a discretely consistent manner
#                 used for to check that the stage cache object is compatible
abstract type AbstractMethod{X, MODE, NS} end

isadjoint(::AbstractMethod{X, MODE}) where {X, MODE} = isadjoint(MODE)
nstages(::AbstractMethod{X, MODE, NS}) where {X, MODE, NS} = NS
mode(::AbstractMethod{X, MODE}) where {X, MODE} = MODE
