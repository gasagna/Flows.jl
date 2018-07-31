
# Four-register implementation of [3R] IMEXRK schemes from section 1.2.3 of CB 2015
function step!(scheme::Scheme{X, :_4R3R},
               sys::System,
               t::Real,
               Δt::Real,
               x::X) where {X}

    # hoist temporaries out
    y   = scheme.cache[1]
    zᴵ  = scheme.cache[2]
    zᴱ  = scheme.cache[3]
    w   = scheme.cache[4]
    tab = scheme.tableau

    # loop over stages
    for k = 1:nstages(tab)
        if k == 1
             y .= x
            zᴵ .= x
            zᴱ .= x
        else
            zᴱ .= y .+ tab[:aᴱ, k, k-1].*Δt.*zᴱ
            if k < nstages(tab)
                y .= x .+ (tab[:aᴵ, k+1, k-1] .- tab[:bᴵ, k-1]).*Δt.*zᴵ .+
                          (tab[:aᴱ, k+1, k-1] .- tab[:bᴱ, k-1])./tab[:aᴱ, k, k-1].*(zᴱ .- y)
            end
            zᴱ .= zᴱ .+ tab[:aᴵ, k, k-1]*Δt.*zᴵ
        end
        A_mul_B!(w, sys, zᴱ)                 # compute w = A*zᴱ then
        ImcA!(sys, tab[:aᴵ, k, k]*Δt, w, zᴵ) # compute z = (I-cA)⁻¹*w in place
        w .= zᴱ .+ tab[:aᴵ, k, k]*Δt.*zᴵ     # w is the temp input, output is zᴱ
        sys(t + tab[:cᴱ, k]*Δt, w, zᴱ)
        x .= x .+ tab[:bᴵ, k].*Δt.*zᴵ .+ tab[:bᴱ, k].*Δt.*zᴱ
    end
    return nothing
end