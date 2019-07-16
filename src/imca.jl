import LinearAlgebra

export ImcA!

"""
    ImcA!(A, c::Real, y, z)

Return `z` that solves the linear problem `(I - c*A)*z = y`, where `c` is a scalar, 
`A` a linear operator and `I` is the identity. 

To use the IMEX schemes in this package, user should add methods to this function for 
their custom types, using an efficient implementation. A default implementation when 
`A` is of type `LinearAlgebra.Diagonal` and `y` and `z` is provided by this package.

# Notes

The name of this function should be read "I-minus-see-A".
"""
ImcA!(A, c::Real, y, z) = 
    error("ImcA! missing implementation for operator `A` of type $(typeof(A))")

"""
    ImcA_mul!(A, c::Real, y, z)

Calculate the product 'z = (I - c*A)*y', where `c` is a scalar, `A` is a linear operator 
and `I` is the identity. 
"""
ImcA_mul!(A, c::Real, y, z) = (mul!(z, A, y); z .= y .- c.*z; z) 


# Provide interface for systems where the stiff operator is
# defined as a `Diagonal` matrix object from Julia Base.
ImcA!(A::LinearAlgebra.Diagonal, c::Real, y::V, z::V) where {V<:AbstractVector} = 
    (z .= y./(1 .- c.*A.diag); z)

LinearAlgebra.mul!(out::V, A::LinearAlgebra.Diagonal, in::V) where {V<:AbstractVector} = 
    (out .= A.diag.*in; out)