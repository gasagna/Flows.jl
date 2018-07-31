export RK4, RK4_Tan, RK4_Adj

# ---------------------------------------------------------------------------- #
# Classical fourth order Runge-Kutta
struct RK4{X, ISCACHE} <: AbstractMethod
    store::NTuple{5, X}
end

RK4(x::X, store::Bool=false) where {X} = RK4{X, store}(ntuple(i->similar(x), 5))

function step!(method::RK4{X, ISCACHE},
                  sys::System,
                    t::Real,
                   Δt::Real,
                    x::X,
              stcache::C) where {X, ISCACHE, C<:AbstractStageCache{4, X}}
    # aliases
    k1, k2, k3, k4, y = method.store

    # stages
    y .= x;              sys(t,        y, k1); ISCACHE && (s1 = copy(y))
    y .= x .+ Δt.*k1./2; sys(t + Δt/2, y, k2); ISCACHE && (s2 = copy(y))
    y .= x .+ Δt.*k2./2; sys(t + Δt/2, y, k3); ISCACHE && (s3 = copy(y))
    y .= x .+ Δt.*k3   ; sys(t + Δt,   y, k4); ISCACHE && (s4 = copy(y))

    # store stages if requested
    ISCACHE && push!(stcache, t, Δt, (s1, s2, s3, s4))

    # wrap up
    x .= x .+ Δt./6 .* (k1 .+ 2.0.*k2 .+ 2.0.*k3 .+ k4)

    return nothing
end

# ---------------------------------------------------------------------------- #
# Linearisation of classical fourth order Runge-Kutta
struct RK4_Tan{X} <: AbstractMethod
    store::NTuple{5, X}
end

# constructor
RK4_Tan(x::X) where {X} = RK4_Tan{X}(ntuple(i->similar(x), 5))

# takes x_{n} and overwrites it with x_{n+1}
function step!(method::RK4_Tan{X},   #
               sys::System,          #
               t::Real,              # the time corresponding to x_{n}
               Δt::Real,             # will be positive
               x::X,
               stages::NTuple{4, X}) where {X}
    # aliases
    k1, k2, k3, k4, y = method.store

    # stages
    y .= x               ; sys(t,        stages[1], y, k1)
    y .= x .+ 0.5.*Δt.*k1; sys(t + Δt/2, stages[2], y, k2)
    y .= x .+ 0.5.*Δt.*k2; sys(t + Δt/2, stages[3], y, k3)
    y .= x .+      Δt.*k3; sys(t + Δt,   stages[4], y, k4)

    # wrap up
    x .= x .+ Δt./6 .* (k1 .+ 2.0.*k2 .+ 2.0.*k3 .+ k4)

    return nothing
end

# ---------------------------------------------------------------------------- #
# Adjoint of linearisation of classical fourth order Runge-Kutta
struct RK4_Adj{X} <: AbstractMethod
    store::NTuple{5, X}
end

# constructor
RK4_Adj(x::X) where {X} = RK4_Adj{X}(ntuple(i->similar(x), 5))

# takes x_{n+1} and overwrites it with x_{n}
function step!(method::RK4_Adj{X},   #
               sys::System,          #
               t::Real,              # the time corresponding to x_{n}
               Δt::Real,             # will be positive
               x::X,
               stages::NTuple{4, X}) where {X}
    # aliases
    j1, j2, j3, j4, y = method.store

    # stages
    y .= x               ; sys(t + Δt  , stages[4], y, j4)
    y .= x .+ 0.5.*Δt.*j4; sys(t + Δt/2, stages[3], y, j3)
    y .= x .+ 0.5.*Δt.*j3; sys(t + Δt/2, stages[2], y, j2)
    y .= x .+      Δt.*j2; sys(t       , stages[1], y, j1)
    x .= x .+ Δt./6 .* (j1 .+ 2.0.*j2 .+ 2.0.*j3 .+ j4)

    return nothing
end