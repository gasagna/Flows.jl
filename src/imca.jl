export ImcA!

"""
    ImcA!(A, c, y, z)

Helper function for solving linear systems associated to 3R2R/4R3R 
schemes. Returns `z` such that `(I-cA)z = y`. For custom types, users
should define a custom method for this function.
"""
function ImcA!(A, c::Real, y::T, z::T) where T end

# Provide interface for systems where the stiff operator is 
# defined as a `Diagonal` matrix object from base Julia. 
ImcA!(A::Diagonal, c::Real, y::T, z::T) where T <: AbstractVector =
    z .= y./(1 .- c.*diag(A))