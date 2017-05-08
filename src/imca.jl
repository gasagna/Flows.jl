export ImcA!

"""
    ImcA!(A, c, y, z)

Helper function for solving linear systems associated to 3R2R/4R3R 
schemes. Returns `z` such that `(I-cA)z = y`. For custom types, users
should define method for this function.
"""
function ImcA!(A, c::Real, y::T, z::T) where T end

# Provide interface for systems where the stiff operator is 
# defined as a `Diagonal` matrix object from base Julia. 
function ImcA!(A::Diagonal, c::Real, y::T, z::T) where T
    length(linearindices(y)) == 
        length(linearindices(z)) == 
            size(A, 1) || throw(DimensionMismatch("wrong input size"))
    d = diag(A) # hoist diagonal
    @simd for i in eachindex(z)
        @inbounds z[i] = y[i]/(1 - c*d[i])
    end
    z
end