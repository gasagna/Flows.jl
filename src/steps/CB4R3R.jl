export CB4R3R4

# type def
struct CB4R3R4{X, NS, TAG, ISADJ} <: AbstractMethod{X, NS, TAG, ISADJ}
    store::NTuple{5, X}
end

# outer constructor
CB4R3R4(x::X, tag::Symbol) where {X} =
    CB4R3R4{X, 6, _tag_map(tag)...}(ntuple(i->similar(x), 5))

# required to cope with buggy julia deepcopy implementation
function Base.deepcopy_internal(x::CB4R3R4, dict::IdDict)
    if !( haskey(dict, x) )
        dict[x] = CB4R3R4(x.store[1], typetag(x))
    end
    return dict[x]
end

# Four-register implementation of [3R] IMEXRK schemes from section 1.2.3 of CB 2015
function step!(method::CB4R3R4{X, NS, :NORMAL},
                  sys::System,
                    t::Real,
                   Δt::Real,
                    x::X,
                    c::C) where {X, NS, C<:Union{Nothing, AbstractStageCache{NS, X}}}

    # hoist temporaries out
    y, zᴵ, zᴱ, w, u   = method.store
    tab = CB4

    # temporary vector for storing stages
    _iscache(C) && (stages = sizehint!(X[], NS))

    # loop over stages
    @inbounds for k = 1:NS
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
        mul!(w, sys, zᴱ)                 # compute w = A*zᴱ then
        ImcA!(sys, tab[:aᴵ, k, k]*Δt, w, zᴵ) # compute z = (I-cA)⁻¹*w in place
        w .= zᴱ .+ tab[:aᴵ, k, k]*Δt.*zᴵ     # w is the temp input, output is zᴱ
        sys(t + tab[:cᴱ, k]*Δt, w, zᴱ); _iscache(C) && (push!(stages, copy(w)))
        x .= x .+ tab[:bᴵ, k].*Δt.*zᴵ .+ tab[:bᴱ, k].*Δt.*zᴱ
    end

    # cache stages if requested
    _iscache(C) && push!(c, t, Δt, tuple(stages...))

    return nothing
end

# Four-register implementation of [3R] IMEXRK schemes from section 1.2.3 of CB 2015
function step!(method::CB4R3R4{X, NS, :LIN, ISADJ},
                  sys::System,
                    t::Real,
                   Δt::Real,
                    x::X,
                store::AbstractStorage) where {X, NS, ISADJ}

    # hoist temporaries out
    y, zᴵ, zᴱ, w, u   = method.store
    tab = CB4

    # modifier for the location of the interpolation
    _m_ = ISADJ == true ? -1 : 1

    # loop over stages
    @inbounds for k = 1:NS
        if k == 1
             y .= x
            zᴵ .= x
            zᴱ .= x
        else
            zᴱ .= y .+ tab[:aᴱ, k, k-1].*_m_.*Δt.*zᴱ
            if k < nstages(tab)
                y .= x .+ (tab[:aᴵ, k+1, k-1] .- tab[:bᴵ, k-1]).*_m_.*Δt.*zᴵ .+
                           (tab[:aᴱ, k+1, k-1] .- tab[:bᴱ, k-1])./tab[:aᴱ, k, k-1].*(zᴱ .- y)
            end
            zᴱ .= zᴱ .+ tab[:aᴵ, k, k-1].*_m_.*Δt.*zᴵ
        end
        mul!(w, sys, zᴱ)                 # compute w = A*zᴱ then
        ImcA!(sys, _m_*tab[:aᴵ, k, k]*Δt, w, zᴵ) # compute z = (I-cA)⁻¹*w in place
        w .= zᴱ .+ tab[:aᴵ, k, k].*_m_.*Δt.*zᴵ     # w is the temp input, output is zᴱ
        sys(t + tab[:cᴱ, k]*Δt, store(u, t + tab[:cᴱ, k]*Δt), w, zᴱ)
        x .= x .+ tab[:bᴵ, k].*_m_.*Δt.*zᴵ .+ tab[:bᴱ, k].*_m_.*Δt.*zᴱ
    end

    return nothing
end