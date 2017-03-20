export IMEXRKScheme, IMEXRK3R2R


# Type definition
# E : Bool, whether we have to calculate the embedded scheme too
# T : IMEXTableau
# S : AbstractVector, the type of the storage
struct IMEXRKScheme{T<:IMEXTableau, S<:AbstractVector, E, CODE}
    storage::Vector{S}
    # internal constructor allocates the storage
    function IMEXRKScheme{T, S, E, CODE}(x::S, N::Int) where {T, S, E, CODE}
        new(S[similar(x) for i = 1:N])
    end
end

# External constructor. Adds one storage for embedded schemes
IMEXRKScheme(tab::IMEXTableau, x::AbstractVector, N::Integer, code::Symbol, embed::Bool=false) =
    IMEXRKScheme{typeof(tab), typeof(x), embed, code}(x, embed ? N+1 : N)

# get the tableau
tableau{T, S, E, CODE}(::Type{IMEXRKScheme{T, S, E, CODE}}) = T

# whether the embedded scheme is active or not
isembedded{T, S, E, CODE}(::Type{IMEXRKScheme{T, S, E, CODE}}) = E


# ~~ IMEXRK3R2R ~~
# Three-register implementation of [2R] IMEXRK schemes from section 1.2.1 of CB 2015
IMEXRK3R2R{T, S, E} = IMEXRKScheme{T, S, E, :_3R2R}
IMEXRK3R2R(tab::IMEXTableau, x::AbstractVector, embed::Bool=false) = 
    IMEXRKScheme(tab, x, 3, :_3R2R, embed)

function _step!{T<:IMEXRK3R2R}(I::Type{T}, g, A, t, Δt, x)
    # extract tableau
    tab = tableau(T)

    # start building expression
    expr_all = Expr(:block)

    # hoist temporaries out
    push!(expr_all.args, :(y = I.storage[1]))
    push!(expr_all.args, :(z = I.storage[2]))
    push!(expr_all.args, :(w = I.storage[3]))

    if isembedded(I)
        push!(expr_all.args, :(x̂ = I.storage[4]))
        push!(expr_all.args, :(x̂ .= x))
    end

    # loop over stages
    for k = 1:nstages(tab)
        expr = Expr(:block)
        aᴵkk   = tab[Val{:aᴵ}, k, k]
        bᴱk    = tab[Val{:bᴱ}, k]
        bᴵk    = tab[Val{:bᴵ}, k]
        cᴱk    = tab[Val{:cᴱ}, k]
        b̂ᴱk    = tab[Val{:b̂ᴱ}, k]
        b̂ᴵk    = tab[Val{:b̂ᴵ}, k]
        if k == 1
            push!(expr.args, :(y .= x))
        else
            aᴵkkm1 = tab[Val{:aᴵ}, k, k-1]
            aᴱkkm1 = tab[Val{:aᴱ}, k, k-1]
            bᴵkm1  = tab[Val{:bᴵ}, k-1]
            bᴱkm1  = tab[Val{:bᴱ}, k-1]
            (aᴵkkm1 - bᴵkm1) == 0 && (aᴱkkm1 - bᴱkm1) == 0 &&
                push!(expr.args, :(@inbounds @simd for i in eachindex(y); y[i] = x[i]; end))
            (aᴵkkm1 - bᴵkm1) == 0 && (aᴱkkm1 - bᴱkm1) != 0 &&
                push!(expr.args, :(@inbounds @simd for i in eachindex(y); y[i] = x[i] + $(aᴱkkm1 - bᴱkm1)*Δt*y[i]; end))
            (aᴵkkm1 - bᴵkm1) != 0 && (aᴱkkm1 - bᴱkm1) == 0 &&
                push!(expr.args, :(@inbounds @simd for i in eachindex(y); y[i] = x[i] + $(aᴵkkm1 - bᴵkm1)*Δt*z[i]; end))
            (aᴵkkm1 - bᴵkm1) != 0 && (aᴱkkm1 - bᴱkm1) != 0 &&
                push!(expr.args, :(@inbounds @simd for i in eachindex(y); y[i] = x[i] + $(aᴵkkm1 - bᴵkm1)*Δt*z[i] + $(aᴱkkm1 - bᴱkm1)*Δt*y[i]; end))
        end
        # compute z = A*y then
        # compute z = (I-cA)⁻¹*(A*y) in place
        push!(expr.args, :(A_mul_B!(A, y, z)))
        push!(expr.args, :(ImcA!(A, $aᴵkk*Δt, z, z)))

        # w is the temporary input @inbounds @simd for g - output in y
        push!(expr.args, :(@inbounds @simd for i in eachindex(w); w[i] = y[i]; end))
        aᴵkk != 0 && push!(expr.args, :(@inbounds @simd for i in eachindex(w); w[i] += $aᴵkk*Δt*z[i]; end))
        push!(expr.args, :(g(t + $cᴱk*Δt, w, y)))

        bᴵk != 0 && push!(expr.args, :(@inbounds @simd for i in eachindex(x); x[i] += $bᴵk*Δt*z[i]; end))
        bᴱk != 0 && push!(expr.args, :(@inbounds @simd for i in eachindex(x); x[i] += $bᴱk*Δt*y[i]; end))

        if isembedded(I)
            b̂ᴵk != 0 && push!(expr.args, :(@inbounds @simd for i in eachindex(x̂); x̂[i] += $b̂ᴵk*Δt*z[i]; end))
            b̂ᴱk != 0 && push!(expr.args, :(@inbounds @simd for i in eachindex(x̂); x̂[i] += $b̂ᴱk*Δt*y[i]; end))
        end
        push!(expr_all.args, expr)
    end
    expr_all
end

# ~~

# each IMEXRKScheme implements a _step! method that will be called here by dispatch
@generated function step!(I::IMEXRKScheme, g, A::AbstractMatrix, t::Real, Δt::Real, x::AbstractVector)
    _step!(I, g, A, t, Δt, x)
end