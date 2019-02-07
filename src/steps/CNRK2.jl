export CNRK2

# ---------------------------------------------------------------------------- #
# Crank-Nicolson/Heun used in Chandler and Kerswell 2013
struct CNRK2{X, TAG, ISADJ} <: AbstractMethod{X, 2, TAG, ISADJ}
    store::NTuple{5, X}
end

# outer constructor
CNRK2(x::X, tag::Symbol) where {X} =
    CNRK2{X, _tag_map(tag)...}(ntuple(i->similar(x), 5))

# required to cope with nuggy julia deepcopy implementation
function Base.deepcopy_internal(x::CNRK2, dict::IdDict)
    if !( haskey(dict, x) )
        dict[x] = CNRK2(x.store[1], typetag(x))
    end
    return dict[x]
end

# ---------------------------------------------------------------------------- #
# Normal time stepping with optional stage caching
function step!(method::CNRK2{X, :NORMAL},
                  sys::System,
                    t::Real,
                   Δt::Real,
                    x::X,
                    c::C) where {X, C<:Union{Nothing, AbstractStageCache{2, X}}}
    # aliases
    k1, k2, k3, k4, k5 = method.store

    # predictor step
    ImcA_mul!(sys, -0.5*Δt, x, k1)
    sys(t, x, k2); _iscache(C) && (s1 = copy(x))
    @all k3 .= k1 .+ Δt.*k2
    ImcA!(sys, 0.5*Δt, k3, k4)

    # corrector
    sys(t, k4, k5); _iscache(C) && (s2 = copy(k4))
    @all k3 .= k1 .+ 0.5.*Δt.*(k2 .+ k5)
    ImcA!(sys, 0.5*Δt, k3, x)

    # store stages if requested
    _iscache(C) && push!(c, t, Δt, (s1, s2))

    return nothing
end