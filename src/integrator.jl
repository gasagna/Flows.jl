export integrator,
       set_Δt!

mutable struct Integrator{S<:System, M<:Scheme, TS<:AbstractTimeStepping}
      system::S  # the system to be integrated
      scheme::M  # the scheme, with storage and RK implementation
    stepping::TS # the time step. Will be replaced by integration options
    function Integrator(system::S,
                        scheme::M,
                        stepping::TS) where {S, M, TS}
        new{S, M, TS}(system, scheme, stepping)
    end
end

# Outer constructor
integrator(g, A, q, scheme::Scheme, stepping::AbstractTimeStepping) = 
    Integrator(System(g, A,             q), scheme, stepping)
integrator(g, A,    scheme::Scheme, stepping::AbstractTimeStepping) = 
    Integrator(System(g, A,       nothing), scheme, stepping)
integrator(g,       scheme::Scheme, stepping::AbstractTimeStepping) = 
    Integrator(System(g, nothing, nothing), scheme, stepping)
# ---------------------------------------------------------------------------- #
# PROPAGATION FUNCTIONS

# ---------------------------------------------------------------------------- #
# PROPAGATION FUNCTION FOR CONSTANT TIME STEPPING
function _propagate!(method::AbstractMethod{X},
                   stepping::TimeStepConstant
                     system::System,
                       span::NTuple{2, Real},
                          z::X,
                        mon::M,
                      cache::C) where {X, 
                                       M<:Union{Void, AbstractMonitor}, 
                                       C<:Union{Void, AbstractCache}}
    # only allow forward integration
    span[1] == span[2] && return z
    span[1]  > span[2] && throw(ArgumentError("invalid integration time span"))

    # define integration times
    ts = LossLessRange(span[1], span[2], stepping.Δt)

    # push initial state monitor
    _ismonitor(M) && push!(mon, ts[1], _state(z))

    # start integration
    for j = 2:length(ts)

        # advance
        step!(method, system, ts[j-1], ts[j]-ts[j-1], z, cache)

        # store solution into monitor
        _ismonitor(M) && push!(mon, ts[j], _state(z))
    end

    return z
end
                                system::System,
                                span::NTuple{2, Real},
                                stepping::TimeStepConstant,
                                z::S,
                                mon::Vararg{AbstractMonitor, N}) where {S, N}
    quote
        $(Expr(:meta, :inline))

        # define integration times
        ts = LossLessRange(span[1], span[2], 
                           span[1] < span[2] ? stepping.Δt : -stepping.Δt)

        # store initial state in monitors
        Base.Cartesian.@nexprs $N i->push!(mon[i], ts[1], _state(z))

        # start integration
        for j = 2:length(ts)

            # advance
            step!(scheme, system, ts[j-1], ts[j]-ts[j-1], z)

            # store solution into monitor
            Base.Cartesian.@nexprs $N i->push!(mon[i], ts[j], _state(z))
        end
        z
    end
end

@generated function _propagate!(scheme::Scheme{S},
                                system::System,
                                span::NTuple{2, Real},
                                hook::AbstractTimeStepFromHook,
                                z::S,
                                mon::Vararg{AbstractMonitor, N}) where {S, N}
    quote
        $(Expr(:meta, :inline))

        # init and final times
        t, T = span

        # store initial state in monitors
        Base.Cartesian.@nexprs $N i->push!(mon[i], t, _state(z))

        # run until condition
        while t != T
            # obtain time step from hook and what the next time is
            t_next, Δt = _next_Δt(t, T, hook(system.g, system.A, z))

            # advance
            step!(scheme, system, t, Δt, z)

            # update
            t = t_next

            # store solution into monitor
            Base.Cartesian.@nexprs $N i->push!(mon[i], t, _state(z))
        end
        z
    end
end

# return next time and time step that allows hitting the end point exactly
function _next_Δt(t, T, Δt::S) where {S<:Real}
    @assert Δt > 0 "negative time step encountered"
    t_next = ifelse(t ≤ T, min(t+Δt, T), max(t-Δt, T))
    return t_next, S(t_next - t)
end