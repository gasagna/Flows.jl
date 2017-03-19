export DiagonalIMEXSystem

# Dynamical system in which the stiff part is a diagonal linear operator  
struct DiagonalIMEXSystem{G, D} <: AbstractIMEXSystem{G, Diagonal}
    g::G
    A::D
end

# The following four functions are required to make DiagonalIMEXSystem
# work with the rest of the code. In general other subtypes of AbstractIMEXSystem
# should define these four functions.

# Return stiff and non stiff parts
stiff(d::DiagonalIMEXSystem) = d.A
nonstiff(d::DiagonalIMEXSystem) = d.g

# Calculate matrix vector product y = A*x
function Base.A_mul_B!(A::Diagonal, x::AbstractVector, y::AbstractVector)
    length(y) == length(x) == size(A, 1) || throw(DimensionMismatch("wrong input size"))
    d = diag(A) # hoist diagonal
    @simd for i in eachindex(x)
        @inbounds y[i] = d[i]*x[i]
    end
    y
end

# Helper function for solving linear systems associated to 3R2R/2R2R schemes
# Returns z such that (I-cA)z = y.
function ImcA!(A::Diagonal, c::Real, y::AbstractVector, z::AbstractVector)
    length(y) == length(z) == size(A, 1) || throw(DimensionMismatch("wrong input size"))
    d = diag(A) # hoist diagonal
    oneA = one(one(eltype(d))*c)
    @simd for i in eachindex(z)
        @inbounds z[i] = y[i]/(oneA - c*d[i])
    end
    z
end