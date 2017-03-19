export AbstractIMEXSystem, stiff, nonstiff, ImcA!

# Abstract type for all dynamical system that have a stiff
# and non stiff component, where the stiff part is a linear operator
# Concrete implementations of this type can be defined for example
# to allocate specific temporaries, if required. Concrete implementation
# must implement:
#
#     A = stiff(::AbstractIMEXSystem)
#     g = nostiff(::AbstractIMEXSystem)
#
# where A will be a subtype of matrix and implements the following methods
#
#    Base.A_mul_B!(A, x::AbstractVector, y::AbstractVector)
#   ImcA!(A, c::Real, y::AbstractVector, z::AbstractVector)
#
# whereas g must be callable, with signature
#    g(t, x, xdot)
abstract type AbstractIMEXSystem{G, M<:AbstractMatrix} end

# concrete subtypes will satisfy this interface
function stiff(::AbstractIMEXSystem) end
function nonstiff(::AbstractIMEXSystem) end
function ImcA!(A::AbstractMatrix, c::Real, y::AbstractVector, z::AbstractVector) end