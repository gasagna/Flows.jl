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
# A must be diagonal
function ImcA!(A::Diagonal, c::Real, y::AbstractVector, z::AbstractVector)
    length(y) == length(z) == size(A, 1) || throw(DimensionMismatch("wrong input size"))
    d = diag(A) # hoist diagonal
    oneA = one(one(eltype(d))*c)
    @simd for i in eachindex(z)
        @inbounds z[i] = y[i]/(oneA - c*d[i])
    end
    z
end

"""
Three-register implementation of [2R] IMEXRK schemes from section 1.2.1 of CB 2015
 
    Arguments
    ---------
    g   : a in-place function with signature g(t, y, ẏ) - explicit part of the rhs
    A   : a diagonal matrix of the implicit part of the rhs
    x   : state vector at time t, will contain solution of main schem at t+Δt
    x̂   : state vector at time t, will contain solution of embedded scheme at t+Δt
    y   : temporary
    z   : temporary
    w   : temporary
    tab : IMEXRK scheme
    t   : current time
    Δt  : time step
"""
@generated function _3R2R(g, A::Diagonal, x::Vector, x̂::Vector, y::Vector, z::Vector,
                             w::Vector, tab::IMEXTableau, t::Real, Δt::Real)
    _3R2R_impl(g, A, x, x̂, y, z, w, tab, t, Δt)
end

function _3R2R_impl(g, A, x, x̂, y, z, w, tab, t, Δt)
    expr_all = Expr(:block)
    # copy x into x̂
    push!(expr_all.args, :(x̂ .= x))
    # convert tableau coefficients to type of storage
    T = eltype(x) 
    for k = 1:nstages(tab)
        expr = Expr(:block)
        aᴵkk   = T(tab[Val{:aᴵ}, k, k])
        bᴱk    = T(tab[Val{:bᴱ}, k])
        b̂ᴱk    = T(tab[Val{:b̂ᴱ}, k])
        bᴵk    = T(tab[Val{:bᴵ}, k])
        b̂ᴵk    = T(tab[Val{:b̂ᴵ}, k]) 
        cᴱk    = T(tab[Val{:cᴱ}, k])
        if k == 1
            push!(expr.args, :(y .= x))
        else
            aᴵkkm1 = T(tab[Val{:aᴵ}, k, k-1])
            aᴱkkm1 = T(tab[Val{:aᴱ}, k, k-1])
            bᴵkm1  = T(tab[Val{:bᴵ}, k-1])
            bᴱkm1  = T(tab[Val{:bᴱ}, k-1])
            (aᴵkkm1 - bᴵkm1) == 0 && (aᴱkkm1 - bᴱkm1) == 0 && 
                push!(expr.args, :(y .= x))
            (aᴵkkm1 - bᴵkm1) == 0 && (aᴱkkm1 - bᴱkm1) != 0 && 
                push!(expr.args, :(y .= x .+ $(aᴱkkm1 - bᴱkm1)*Δt.*y))
            (aᴵkkm1 - bᴵkm1) != 0 && (aᴱkkm1 - bᴱkm1) == 0 && 
                push!(expr.args, :(y .= x .+ $(aᴵkkm1 - bᴵkm1)*Δt.*z))
            (aᴵkkm1 - bᴵkm1) != 0 && (aᴱkkm1 - bᴱkm1) != 0 && 
                push!(expr.args, :(y .= x .+ $(aᴵkkm1 - bᴵkm1)*Δt.*z .+ $(aᴱkkm1 - bᴱkm1)*Δt.*y))
        end
        # compute z = A*y then
        # compute z = (I-cA)⁻¹*(A*y) in place
        push!(expr.args, :(A_mul_B!(A, y, z)))
        push!(expr.args, :(ImcA!(A, $aᴵkk*Δt, z, z)))

        # w is the temporary input for g - output in y
        push!(expr.args, :(w .= y))
        aᴵkk != 0 && push!(expr.args, :(w .+= ($aᴵkk*Δt).*z))
        push!(expr.args, :(g(t + $cᴱk*Δt, w, y)))

        bᴵk != 0 && push!(expr.args, :(x .+= ($bᴵk*Δt).*z))
        bᴱk != 0 && push!(expr.args, :(x .+= ($bᴱk*Δt).*y))

        b̂ᴵk != 0 && push!(expr.args, :(x̂ .+= ($b̂ᴵk*Δt).*z))
        b̂ᴱk != 0 && push!(expr.args, :(x̂ .+= ($b̂ᴱk*Δt).*y))
        push!(expr_all.args, expr)
    end
    expr_all
end