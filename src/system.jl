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
# where A implements the following methods
#
#    Base.A_mul_B!(A, x::AbstractVector, y::AbstractVector)
#   ImcA!(A, c::Real, y::AbstractVector, z::AbstractVector)
#
# whereas g must be callable, with signature
#    g(t, x, xdot)
abstract type AbstractIMEXSystem{G, A} end

# Concrete subtypes will satisfy this interface

"""
    stiff(f)

Return the stiff part of the system.
"""
function stiff(::AbstractIMEXSystem) end

"""
    nonstiff(f)

Return the nonstiff part of the system.
"""    
function nonstiff(::AbstractIMEXSystem) end

"""
    ImcA!(A, c, y, z)

Helper function for solving linear systems associated to 3R2R/4R3R 
schemes. Returns `z` such that `(I-cA)z = y`. 
"""
function ImcA!(A, c::Real, y::AbstractVector, z::AbstractVector) end