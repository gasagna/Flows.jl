# Type parameters:
# - X           : the state type
# - ISADJOINT   : true to denote adjoint integration
# - NS          : number of internal stages that are required to be saved to then 
#                 integrate linearised equations in a discretely consistent manner
#                 used for to check that the stage cache object is compatible
abstract type AbstractMethod{X, ISADJOINT, NS} end

#
isadjoint(::AbstractMethod{X, ISADJOINT}) where {X, ISADJOINT} = ISADJOINT

# number of internal stages that are required to be saved to then 
# integrate linearised equations in a discretely consistent manner
nstages(::AbstractMethod{X, ISADJOINT, NS}) where {X, ISADJOINT, NS} = NS