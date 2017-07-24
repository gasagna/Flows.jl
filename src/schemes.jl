export IMEXRKScheme, IMEXRK3R2R, IMEXRK4R3R, isembedded, step

"""
    Type that holds information about the particular scheme and 
    implementation, manages the additional storage required for the
    scheme and for the lower-order embedded scheme if that is available.

    One should not use this type directly. The aliases, providing a 
    specific implementation should be used instead.

    Type parameters
    ---------------
    T : IMEXTableau - the tableau of the scheme
    S : the type of the storage
    E : Bool, whether we have to calculate the embedded scheme too. If false
              code blocks that calculate the next state based on the embedded
              scheme are elided from generation.
    CODE : Symbol, the code of the particular implementation used, e.g. _3R2R
              for the three register implementation of 2R scheme. This is
              the only parameter used for dispatch, because for each 
              implementation, i.e. for each CODE there is a particular 
              definition of the method _step!.

"""
struct IMEXRKScheme{T<:IMEXTableau, S, E, CODE}
    storage::Vector{S}
    # internal constructor allocates the storage
    function IMEXRKScheme{T, S, E, CODE}(N::Int, x::S) where {T, S, E, CODE}
        new(S[similar(x) for i = 1:N])
    end
end

"""
    IMEXRKScheme(tab, N, code, embed, x, q...)

Construct an IMEX scheme, using tableau `tab`, where `N` items similar to `x` 
are used for storage. The scheme is parametrised by `code`, that is used for 
dispatching to a particular implementation of the `_step!` method. If `embed`
is true, the embedded scheme is activated, and the next state calculated by the
embedded scheme can be found in the last storage element of the scheme, after
`_step!` has been called to go the the next state.
"""
function IMEXRKScheme(tab::IMEXTableau, N::Integer, code::Symbol, embed::Bool, x, q...)
    z = aug_state(x, q...)
    IMEXRKScheme{typeof(tab), typeof(z), embed, code}(embed ? N+1 : N, z)
end

""" Get the tableau of the scheme """
tableau{T, S, E, CODE}(::Type{IMEXRKScheme{T, S, E, CODE}}) = T

""" Whether the embedded scheme is active or not """
isembedded{T, S, E, CODE}(::Type{IMEXRKScheme{T, S, E, CODE}}) = E
isembedded{T, S, E, CODE}(::IMEXRKScheme{T, S, E, CODE}) = E


# ~~ IMEXRK3R2R ~~
# Three-register implementation of [2R] IMEXRK schemes from section 1.2.1 of CB 2015
IMEXRK3R2R{T, S, E} = IMEXRKScheme{T, S, E, :_3R2R}
IMEXRK3R2R(tab::IMEXTableau, embed::Bool, x, q...) = 
    IMEXRKScheme(tab, 3, :_3R2R, embed, x, q...)

function _step!{T<:IMEXRK3R2R}(I::Type{T}, g, A, t, Δt, x)
    # extract tableau
    tab = tableau(T)

    # start building expression
    expr_all = Expr(:block)

    # hoist temporaries out
    push!(expr_all.args, :(y = I.storage[1]))
    push!(expr_all.args, :(z = I.storage[2]))
    push!(expr_all.args, :(w = I.storage[3]))

    # code will be generated for state and quadrature
    isQuad = x <: AugmentedState

    if isembedded(I)
        push!(expr_all.args, :(x̂ = I.storage[4]))
        push!(expr_all.args, broadcast2fields(Val{isQuad}, :(@over_i x̂[i] = x[i])))
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
            push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i y[i] = x[i])))
        else
            aᴵkkm1 = tab[Val{:aᴵ}, k, k-1]
            aᴱkkm1 = tab[Val{:aᴱ}, k, k-1]
            bᴵkm1  = tab[Val{:bᴵ}, k-1]
            bᴱkm1  = tab[Val{:bᴱ}, k-1]
            (aᴵkkm1 - bᴵkm1) == 0 && (aᴱkkm1 - bᴱkm1) == 0 &&
                push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i y[i] = x[i])))
            (aᴵkkm1 - bᴵkm1) == 0 && (aᴱkkm1 - bᴱkm1) != 0 &&
                push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i y[i] = x[i] + $(aᴱkkm1 - bᴱkm1)*Δt*y[i])))
            (aᴵkkm1 - bᴵkm1) != 0 && (aᴱkkm1 - bᴱkm1) == 0 &&
                push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i y[i] = x[i] + $(aᴵkkm1 - bᴵkm1)*Δt*z[i])))
            (aᴵkkm1 - bᴵkm1) != 0 && (aᴱkkm1 - bᴱkm1) != 0 &&
                push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i y[i] = x[i] + $(aᴵkkm1 - bᴵkm1)*Δt*z[i] + $(aᴱkkm1 - bᴱkm1)*Δt*y[i])))
        end
        # compute z = A*y then
        # compute z = (I-cA)⁻¹*(A*y) in place
        push!(expr.args, :(A_mul_B!(z, A, y)))
        push!(expr.args, :(ImcA!(A, $aᴵkk*Δt, z, z)))

        # w is the temporary input for g - output in y
        push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i w[i] = y[i])))
        aᴵkk != 0 && push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i w[i] += $aᴵkk*Δt*z[i])))
        push!(expr.args, :(g(t + $cᴱk*Δt, w, y)))

        bᴵk != 0 && push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i x[i] += $bᴵk*Δt*z[i])))
        bᴱk != 0 && push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i x[i] += $bᴱk*Δt*y[i])))

        if isembedded(I)
            b̂ᴵk != 0 && push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i x̂[i] += $b̂ᴵk*Δt*z[i])))
            b̂ᴱk != 0 && push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i x̂[i] += $b̂ᴱk*Δt*y[i])))
        end
        push!(expr_all.args, expr)
    end
    expr_all
end

# ~~ IMEXRK4R3R ~~
# Four-register implementation of [3R] IMEXRK schemes from section 1.2.3 of CB 2015
IMEXRK4R3R{T, S, E} = IMEXRKScheme{T, S, E, :_4R3R}
IMEXRK4R3R(tab::IMEXTableau, embed::Bool, x, q...) = 
    IMEXRKScheme(tab, 4, :_4R3R, embed, x, q...)

function _step!{T<:IMEXRK4R3R}(I::Type{T}, g, A, t, Δt, x)
    # extract tableau
    tab = tableau(T)

    # start building expression
    expr_all = Expr(:block)

    # hoist temporaries out
    push!(expr_all.args, :(y  = I.storage[1]))
    push!(expr_all.args, :(zᴵ = I.storage[2]))
    push!(expr_all.args, :(zᴱ = I.storage[3]))
    push!(expr_all.args, :(w  = I.storage[4]))

    # code will be generated for state and quadrature
    isQuad = x <: AugmentedState

    if isembedded(I)
        push!(expr_all.args, :(x̂ = I.storage[5]))
        push!(expr_all.args, broadcast2fields(Val{isQuad}, :(@over_i x̂[i] = x[i])))
    end

    # loop over stages
    for k = 1:nstages(tab)
        expr = Expr(:block)
        aᴵkk = tab[Val{:aᴵ}, k, k]
        bᴱk  = tab[Val{:bᴱ}, k]
        bᴵk  = tab[Val{:bᴵ}, k]
        cᴱk  = tab[Val{:cᴱ}, k]
        b̂ᴱk  = tab[Val{:b̂ᴱ}, k]
        b̂ᴵk  = tab[Val{:b̂ᴵ}, k]
        if k == 1
            push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i  y[i] = x[i])))
            push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i zᴵ[i] = x[i])))
            push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i zᴱ[i] = x[i])))
        else
            aᴱkkm1 = tab[Val{:aᴱ}, k, k-1]
            aᴱkkm1 == 0 && push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i zᴱ[i] = y[i])))
            aᴱkkm1 != 0 && push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i zᴱ[i] = y[i] + $aᴱkkm1*Δt*zᴱ[i])))
            if k < nstages(tab)
                aᴵkpkm1 = tab[Val{:aᴵ}, k+1, k-1]
                aᴱkpkm1 = tab[Val{:aᴱ}, k+1, k-1]
                bᴵkm1   = tab[Val{:bᴵ}, k-1]
                bᴱkm1   = tab[Val{:bᴱ}, k-1]
                
                c1 =  aᴵkpkm1 - bᴵkm1
                c2 = (aᴱkpkm1 - bᴱkm1)/aᴱkkm1
                c1 == 0 && c2 == 0 && 
                    push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i y[i] = x[i])))
                c1 == 0 && c2 != 0 && 
                    push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i y[i] = x[i] + $c2*(zᴱ[i] - y[i]))))
                c1 != 0 && c2 == 0 && 
                    push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i y[i] = x[i] + $c1*Δt*zᴵ[i])))
                c1 != 0 && c2 != 0 && 
                    push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i y[i] = x[i] + $c1*Δt*zᴵ[i] + $c2*(zᴱ[i] - y[i]))))
            end
            aᴵkkm1 = tab[Val{:aᴵ}, k, k-1]
            aᴵkkm1 != 0 && push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i zᴱ[i] = zᴱ[i] + $aᴵkkm1*Δt*zᴵ[i])))
        end
        # compute w = A*zᴱ then
        # compute z = (I-cA)⁻¹*w in place
        push!(expr.args, :(A_mul_B!(w, A, zᴱ)))
        push!(expr.args, :(ImcA!(A, $aᴵkk*Δt, w, zᴵ)))

        # w is the temporary input for g - output in zᴱ
        push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i w[i] = zᴱ[i])))
        aᴵkk != 0 && push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i w[i] += $aᴵkk*Δt*zᴵ[i])))

        push!(expr.args, :(g(t + $cᴱk*Δt, w, zᴱ)))

        bᴵk != 0 && push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i x[i] += $bᴵk*Δt*zᴵ[i])))
        bᴱk != 0 && push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i x[i] += $bᴱk*Δt*zᴱ[i])))

        if isembedded(I)
            b̂ᴵk != 0 && push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i x̂[i] += $b̂ᴵk*Δt*zᴵ[i])))
            b̂ᴱk != 0 && push!(expr.args, broadcast2fields(Val{isQuad}, :(@over_i x̂[i] += $b̂ᴱk*Δt*zᴱ[i])))
        end
        push!(expr_all.args, expr)
    end
    expr_all
end

# each IMEXRKScheme implements a _step! method that will be called here by dispatch
@generated function step!{T, S, E, CODE}(I::IMEXRKScheme{T, S, E, CODE}, 
                                         g, 
                                         A, 
                                         t::Real, Δt::Real, x::S,)
    _step!(I, g, A, t, Δt, x)
end