export IMEXRK3R2R, IMEXRK4R3R, isembedded, tableau, nstorage

# Type that holds information about the particular implementation
# of an IMEXRK method. Concrete subtypes are primarily used for 
# code generation and dispatch in the generated step! functions.

# Type parameters
# ---------------
# T     : IMEXTableau - the tableau of the scheme
# Embed : Bool, whether we have to calculate the embedded scheme too. If false
#         code blocks that calculate the next state based on the embedded
#         scheme are elided from generation.
# M     : the number of storage areas to be allocated for the implementation
abstract type AbstractIMEXRKImplementation{Tab, Embed, M} end

# all information is in the type
tableau(   ::Type{<:AbstractIMEXRKImplementation{Tab, Embed, M}}) where {Tab, Embed, M} = Tab
isembedded(::Type{<:AbstractIMEXRKImplementation{Tab, Embed, M}}) where {Tab, Embed, M} = Embed
nstorage(  ::Type{<:AbstractIMEXRKImplementation{Tab, Embed, M}}) where {Tab, Embed, M} = Embed == true ? M+1 : M

# works on instances as well
tableau(   ::T) where T = tableau(T)
isembedded(::T) where T = isembedded(T)
nstorage(  ::T) where T = nstorage(T)

# Define types for specific implementations, for which code will be generated
struct IMEXRK3R2R{Tab, Embed} <: AbstractIMEXRKImplementation{Tab, Embed, 3} end
struct IMEXRK4R3R{Tab, Embed} <: AbstractIMEXRKImplementation{Tab, Embed, 4} end

# Constructors for these implementations
IMEXRK3R2R(Tab::IMEXTableau, embed::Bool) = IMEXRK3R2R{Tab, embed}()
IMEXRK4R3R(Tab::IMEXTableau, embed::Bool) = IMEXRK4R3R{Tab, embed}()

# Three-register implementation of [2R] IMEXRK schemes from section 1.2.1 of CB 2015
function _codegen(::IMEXRK3R2R{Tab, Embed}) where {Tab, Embed}

    # start building expression
    expr_all = Expr(:block)

    # get storage
    push!(expr_all.args, :(store = scheme.storage))

    # hoist temporaries out
    push!(expr_all.args, :(y = store[1]))
    push!(expr_all.args, :(z = store[2]))
    push!(expr_all.args, :(w = store[3]))

    if Embed == true
        push!(expr_all.args, :(xe = store[4]))
        push!(expr_all.args, :(xe .= x))
    end

    # loop over stages
    for k = 1:nstages(Tab)
        expr = Expr(:block)

        aᴵkk   = Tab[:aᴵ, k, k]
        bᴱk    = Tab[:bᴱ, k]
        bᴵk    = Tab[:bᴵ, k]
        cᴱk    = Tab[:cᴱ, k]
        eᴱk    = Tab[:eᴱ, k]
        eᴵk    = Tab[:eᴵ, k]
        if k == 1
            push!(expr.args, :(y .= x))
        else
            aᴵkkm1 = Tab[:aᴵ, k, k-1]
            aᴱkkm1 = Tab[:aᴱ, k, k-1]
            bᴵkm1  = Tab[:bᴵ, k-1]
            bᴱkm1  = Tab[:bᴱ, k-1]
            push!(expr.args, :(y .= x .+ $(aᴵkkm1 - bᴵkm1)*Δt.*z .+ $(aᴱkkm1 - bᴱkm1)*Δt.*y))
        end
        # compute z = A*y then
        # compute z = (I-cA)⁻¹*(A*y) in place
        push!(expr.args, :(A_mul_B!(z, A, y)))
        push!(expr.args, :(ImcA!(A, $aᴵkk*Δt, z, z)))

        # w is the temporary input for g - output in y
        push!(expr.args, :(w .= y .+ $aᴵkk*Δt.*z))
        push!(expr.args, :(g(t + $cᴱk*Δt, w, y)))
        push!(expr.args, :(x .= x .+ $bᴵk*Δt.*z .+ $bᴱk*Δt.*y))

        # add code for embedded implementations
        Embed == true && push!(expr.args, :(xe .= xe .+ $eᴵk*Δt.*z .+ $eᴱk*Δt.*y))

        # add stage to overall expression
        push!(expr_all.args, expr)
    end
    expr_all
end

# Four-register implementation of [3R] IMEXRK schemes from section 1.2.3 of CB 2015
function _codegen(::IMEXRK4R3R{Tab, Embed}) where {Tab, Embed}
    # start building expression
    expr_all = Expr(:block)

    # get storage
    push!(expr_all.args, :(store = scheme.storage))

    # hoist temporaries out
    push!(expr_all.args, :(y  = store[1]))
    push!(expr_all.args, :(zᴵ = store[2]))
    push!(expr_all.args, :(zᴱ = store[3]))
    push!(expr_all.args, :(w  = store[4]))

    if Embed == true
        push!(expr_all.args, :(xe  = store[5]))
        push!(expr_all.args, :(xe .= x))
    end

    # loop over stages
    for k = 1:nstages(Tab)
        expr = Expr(:block)
        aᴵkk = Tab[:aᴵ, k, k]
        bᴱk  = Tab[:bᴱ, k]
        bᴵk  = Tab[:bᴵ, k]
        cᴱk  = Tab[:cᴱ, k]
        eᴱk  = Tab[:eᴱ, k]
        eᴵk  = Tab[:eᴵ, k]
        if k == 1
            push!(expr.args, :( y .= x))
            push!(expr.args, :(zᴵ .= x))
            push!(expr.args, :(zᴱ .= x))
        else
            aᴱkkm1 = Tab[:aᴱ, k, k-1]
            push!(expr.args, :(zᴱ .= y .+ $aᴱkkm1*Δt.*zᴱ))
            if k < nstages(Tab)
                aᴵkpkm1 = Tab[:aᴵ, k+1, k-1]
                aᴱkpkm1 = Tab[:aᴱ, k+1, k-1]
                bᴵkm1   = Tab[:bᴵ, k-1]
                bᴱkm1   = Tab[:bᴱ, k-1]
                c1 =  aᴵkpkm1 - bᴵkm1
                c2 = (aᴱkpkm1 - bᴱkm1)/aᴱkkm1
                push!(expr.args, :(y .= x .+ $c1*Δt.*zᴵ .+ $c2.*(zᴱ .- y)))
            end
            aᴵkkm1 = Tab[:aᴵ, k, k-1]
            push!(expr.args, :(zᴱ .= zᴱ .+ $aᴵkkm1*Δt.*zᴵ))
        end
        # compute w = A*zᴱ then
        # compute z = (I-cA)⁻¹*w in place
        push!(expr.args, :(A_mul_B!(w, A, zᴱ)))
        push!(expr.args, :(ImcA!(A, $aᴵkk*Δt, w, zᴵ)))

        # w is the temporary input for g - output in zᴱ
        push!(expr.args, :(w .= zᴱ .+ $aᴵkk*Δt.*zᴵ))
        push!(expr.args, :(g(t + $cᴱk*Δt, w, zᴱ)))
        push!(expr.args, :(x .= x .+ $bᴵk*Δt.*zᴵ .+ $bᴱk*Δt.*zᴱ))
        
        # add code for embedded implementations
        Embed == true && push!(expr.args, :(xe .= xe .+ $eᴵk*Δt.*zᴵ .+ $eᴱk*Δt.*zᴱ))

        # add stage to overall expression
        push!(expr_all.args, expr)
    end
    expr_all
end