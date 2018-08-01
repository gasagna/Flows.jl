export ImcA!

"""
    ImcA!(A, c, y, z)

Helper function for solving linear systems associated to 3R2R/4R3R
schemes. Returns `z` such that `(I-cA)z = y`. For custom types, users
should define a custom method for this function.
"""
ImcA!(A, c::Real, y, z) = error("missing implementation")

# Provide interface for systems where the stiff operator is
# defined as a `Diagonal` matrix object from Julia Base. Custom
# types are supposed to have defined broadcasting operations
# for the dot notation
ImcA!( A::Diagonal, c::Real, y, z) = z .= y./(1 .- c.*diag(A))
ImcAt!(A::Diagonal, c::Real, y, z) = z .= y./(1 .- c.*diag(A))

Base.A_mul_B!( out, A::Diagonal, in) = out .= diag(A).*in
Base.At_mul_B!(out, A::Diagonal, in) = out .= diag(A).*in