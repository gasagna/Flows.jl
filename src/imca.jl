import LinearAlgebra: mul!, Diagonal

export ImcA!, ImcA_mul!

"""
    ImcA!(A, c, y, z)

Returns `z` such that `(I-cA)z = y`.
"""
ImcA!(A, c::Real, y, z) = error("missing implementation")

"""
    ImcA_mul!(A, c, y, z)

Calculate the matrix vector product 'z = (I - c*A)y', where 'c' is a 
scalar, 'A' and operator, and 'I' the identity operator.
"""
ImcA_mul!(A, c::Real, y, z) = (mul!(z, A, y); @all z .= y .- c.*z; z) 


# Provide interface for systems where the stiff operator is
# defined as a `Diagonal` matrix object from Julia Base.
ImcA!( A::Diagonal, c::Real, y, z) = z .= y./(1 .- c.*A.diag)
mul!( out, A::Diagonal, in) = out .= A.diag.*in