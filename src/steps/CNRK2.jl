export CNRK2

# ---------------------------------------------------------------------------- #
# Crank-Nicolson/Heun used in Chandler and Kerswell 2013
struct CNRK2{X, MODE, ISADJOINT, NX} <: AbstractMethod{X, MODE, ISADJOINT, 2}
    store::NX
    CNRK2{X, MODE, ISADJOINT}(store::NX) where {X, MODE, ISADJOINT, N, NX<:NTuple{N, X}} = 
        new{X, MODE, ISADJOINT, NX}(store)
end

# outer constructor
CNRK2(x::X) where {X} = 
    CNRK2{X, NormalMode, false}(ntuple(i->similar(x), 5))

CNRK2(x::X, ::ContinuousMode, isadjoint::Bool=false) where {X} = 
    CNRK2{X, ContinuousMode, isadjoint}(ntuple(i->similar(x), 5))

CNRK2(x::X, ::DiscreteMode, isadjoint::Bool=false) where {X} = 
    CNRK2{X, DiscreteMode, isadjoint}(ntuple(i->similar(x), 5))

# required to cope with nuggy julia deepcopy implementation
function Base.deepcopy_internal(x::CNRK2,
                             dict::IdDict)
    if !( haskey(dict, x) )
        dict[x] = CNRK2(x.store[1], mode(x), isadjoint(x))
    end
    return dict[x]
end

# ---------------------------------------------------------------------------- #
# Normal time stepping with optional stage caching
function step!(method::CNRK2{X, NormalMode},
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
    k3 .= k1 .+ Δt.*k2
    ImcA!(sys, 0.5*Δt, k3, k4)

    # corrector
    sys(t+Δt, k4, k5); _iscache(C) && (s2 = copy(k4))
    k3 .= k1 .+ 0.5.*Δt.*(k2 .+ k5)
    ImcA!(sys, 0.5*Δt, k3, x)

    # store stages if requested
    _iscache(C) && push!(c, t, Δt, (s1, s2))

    return nothing
end

# ---------------------------------------------------------------------------- #
# Continuous time stepping for linearised/adjoint equations with interpolation
# from an `AbstractStorage` object for the evaluation of the linear operator.
function step!(method::CNRK2{X, ContinuousMode, ISADJOINT},
                  sys::System,
                    t::Real,
                   Δt::Real,
                    x::X,
                store::AbstractStorage) where {X, ISADJOINT}

    # modifier for the location of the interpolation
    _m_ = ISADJOINT == true ? -1 : 1

    # aliases
    k1, k2, k3, k4, k5 = method.store
    ImcA_mul!(sys, -0.5*_m_*Δt, x, k1)
    sys(t, store(k5, t), x, k2)
    k3 .= k1 .+ _m_.*Δt.*k2
    ImcA!(sys, 0.5*_m_*Δt, k3, k4)
    sys(t + Δt, store(k3, t + Δt), k4, k5)
    k3 .= k1 .+ 0.5.*_m_.*Δt.*(k2 .+ k5)
    ImcA!(sys, 0.5*_m_*Δt, k3, x)
    return nothing
end

# ---------------------------------------------------------------------------- #
# Forward linearised method takes x_{n} and overwrites it with x_{n+1}
# Adjoint linearised method takes x_{n+1} and overwrites it with x_{n}
function step!(method::CNRK2{X, DiscreteMode, ISADJOINT},
                  sys::System,
                    t::Real,
                   Δt::Real,
                    x::X,
               stages::NTuple{2, X}) where {X, ISADJOINT}
    # aliases
    k1, k2, k3, k4, k5 = method.store

    if ISADJOINT
        ImcA!(sys, 0.5*Δt, x, k1)
        sys(t + Δt, stages[2], k1, k2)
        k2 .= k2.*Δt./2
        ImcA!(sys, 0.5*Δt, k2, k3)
        k4 .= Δt./2.0.*k1 .+ Δt.*k3
        sys(t, stages[1], k4, k5)
        k2 .= k1 .+ k3
        ImcA_mul!(sys, -0.5*Δt, k2, k3)
        x .= k3 .+ k5
    else
        ImcA_mul!(sys, -0.5*Δt, x, k1)
        sys(t, stages[1], x, k2)
        k3 .= k1 .+ Δt.*k2
        ImcA!(sys, 0.5*Δt, k3, k4)
        sys(t+Δt, stages[2], k4, k5)
        k3 .= k1 .+ 0.5.*Δt.*(k2 .+ k5)
        ImcA!(sys, 0.5*Δt, k3, x)
    end

    return nothing
end
