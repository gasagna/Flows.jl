export ContinuousMode, DiscreteMode

# Type parameters:
# - X           : the state type
# - MODE        : integration mode: norma, linear continuous, lienar discrete
# - ISADJOINT   : true to denote adjoint integration
# - NS          : number of internal stages that are required to be saved to then 
#                 integrate linearised equations in a discretely consistent manner
#                 used for to check that the stage cache object is compatible
abstract type AbstractMethod{X, MODE, ISADJOINT, NS} end

# whether the scheme is for the adjoint integration or not
isadjoint(::AbstractMethod{X, MODE, ISADJOINT}) where {X, MODE, ISADJOINT} = ISADJOINT

# number of internal stages that are required to be saved to then 
# integrate linearised equations in a discretely consistent manner
nstages(::AbstractMethod{X, MODE, ISADJOINT, NS}) where {X, MODE, ISADJOINT, NS} = NS

# number of internal stages that are required to be saved to then 
# integrate linearised equations in a discretely consistent manner
mode(::AbstractMethod{X, MODE, ISADJOINT, NS}) where {X, MODE, ISADJOINT, NS} = MODE()

# tag objects to specify the integration mode
struct NormalMode end
struct ContinuousMode end
struct DiscreteMode end