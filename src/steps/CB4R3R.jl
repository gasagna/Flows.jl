export CB4R3R4

# type def
struct CB4R3R4{X, NS, TAG, ISADJ} <: AbstractMethod{X, NS, TAG, ISADJ}
    store::NTuple{4, X}
end

# outer constructor
CB4R3R4(x::X, tag::Symbol) where {X} =
    CB4R3R4{X, 6, _tag_map(tag)...}(ntuple(i->similar(x), 4))


# Four-register implementation of [3R] IMEXRK schemes from section 1.2.3 of CB 2015
function step!(method::CB4R3R4{X, NS, :NORMAL},
                  sys::System,
                    t::Real,
                   Δt::Real,
                    x::X,
                    c::C) where {X, NS, C<:Union{Nothing, AbstractStageCache{NS, X}}}

    # hoist temporaries out
    y, zᴵ, zᴱ, w   = method.store
    tab = CB4

    # temporary vector for storing stages
    _iscache(C) && (stages = sizehint!(X[], NS))

    # loop over stages
    @inbounds for k = 1:NS
        if k == 1
            @all  y .= x
            @all zᴵ .= x
            @all zᴱ .= x
        else
            @all zᴱ .= y .+ tab[:aᴱ, k, k-1].*Δt.*zᴱ
            if k < nstages(tab)
                @all y .= x .+ (tab[:aᴵ, k+1, k-1] .- tab[:bᴵ, k-1]).*Δt.*zᴵ .+
                               (tab[:aᴱ, k+1, k-1] .- tab[:bᴱ, k-1])./tab[:aᴱ, k, k-1].*(zᴱ .- y)
            end
            @all zᴱ .= zᴱ .+ tab[:aᴵ, k, k-1]*Δt.*zᴵ
        end
        mul!(w, sys, zᴱ)                 # compute w = A*zᴱ then
        ImcA!(sys, tab[:aᴵ, k, k]*Δt, w, zᴵ) # compute z = (I-cA)⁻¹*w in place
        @all w .= zᴱ .+ tab[:aᴵ, k, k]*Δt.*zᴵ     # w is the temp input, output is zᴱ
        sys(t + tab[:cᴱ, k]*Δt, w, zᴱ); _iscache(C) && (push!(stages, copy(w)))
        @all x .= x .+ tab[:bᴵ, k].*Δt.*zᴵ .+ tab[:bᴱ, k].*Δt.*zᴱ
    end

    # cache stages if requested
    _iscache(C) && push!(c, t, Δt, tuple(stages...))

    return nothing
end