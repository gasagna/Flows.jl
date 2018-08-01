export flow

# ---------------------------------------------------------------------------- #
# The Flow type!
struct Flow{S<:System, M<:AbstractMethod, TS<:AbstractTimeStepping}
      sys::S  # the system to be integrated
     meth::M  # the method, with storage, implementation and time stepping
    tstep::TS # the method used for time stepping
    Flow(sys::S, m::M, ts::TS) where {S, M, TS} = new{S, M, TS}(sys, m, ts)
end

# main entry point specifying explicit, implicit and quadrature parts
flow(g, A, q, m::AbstractMethod, ts::AbstractTimeStepping) =
     Flow(System(g, A, q), m, ts)

flow(g, A,    m::AbstractMethod, ts::AbstractTimeStepping) =
     Flow(System(g, A, nothing), m, ts)

flow(g,       m::AbstractMethod, ts::AbstractTimeStepping) =
     Flow(System(g, nothing, nothing), m, ts)

# ---------------------------------------------------------------------------- #
# FLOWS ARE CALLABLE OBJECTS: THIS IS THE MAIN INTERFACE

const _Span = NTuple{2, Real}

@inline (I::Flow)(x, span::_Span) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span), x, nothing, nothing)

@inline (I::Flow)(x, span::_Span, m::AbstractMonitor) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span), x, m, nothing)

@inline (I::Flow)(x, span::_Span, c::AbstractCache) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span), x, nothing, c)

@inline (I::Flow)(x, span::_Span, c::AbstractCache, m::AbstractMonitor) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span), x, m, c)

@inline (I::Flow)(x, span::_Span, m::AbstractMonitor, c::AbstractCache) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span), x, m, c)

# same with quadrature
@inline (I::Flow)(x, q, span::_Span) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span),
                augment(x, q), nothing, nothing)

@inline (I::Flow)(x, q, span::_Span, m::AbstractMonitor) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span),
                augment(x, q), nothing, m)

@inline (I::Flow)(x, q, span::_Span, c::AbstractCache) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span),
                augment(x, q), c, nothing)

@inline (I::Flow)(x, q, span::_Span, c::AbstractCache, m::AbstractMonitor) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span),
                augment(x, q), c, m)

@inline (I::Flow)(x, q, span::_Span, m::AbstractMonitor, c::AbstractCache) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span),
                augment(x, q), c, m)


# ---------------------------------------------------------------------------- #
# PROPAGATION FUNCTIONS & UTILS

# check span and/or return
macro _checkspan(span, z)
    quote
       $(esc(span))[1] ==$(esc(span))[2] && return $z
       $(esc(span))[1]  >$(esc(span))[2] && throw(ArgumentError("invalid span"))
    end
end

# ---------------------------------------------------------------------------- #
# CONSTANT FORWARD TIME STEPPING, ONLY FOR STATE EQUATIONS
function _propagate!(method::AbstractMethod{Z, NS, :NL},
                   stepping::TimeStepConstant,
                     system::System,
                       span::_Span,
                          z::Z,
                        mon::M,
                      cache::C) where {Z, NS,
                                       M<:Union{Void, AbstractMonitor},
                                       C<:Union{Void, AbstractCache}}
    # check span is sane
    @_checkspan(span, z)

    # define integration times
    ts = LossLessRange(span[1], span[2], stepping.Δt)

    # push initial state monitor
    _ismonitor(M) && push!(mon, ts[1], _state(z))

    # start integration
    for j = 2:length(ts)
        step!(method, system, ts[j-1], ts[j]-ts[j-1], z, cache)
        _ismonitor(M) && push!(mon, ts[j], _state(z))
    end

    return z
end

# ---------------------------------------------------------------------------- #
# PROPAGATION BASED ON SYSTEMS HOOK: ONLY FOR STATE EQUATIONS
function _propagate!(method::AbstractMethod{Z, NS, :NL},
                       hook::AbstractTimeStepFromHook,
                     system::System,
                       span::_Span,
                          z::Z,
                        mon::M,
                      cache::C) where {Z, NS,
                                       M<:Union{Void, AbstractMonitor},
                                       C<:Union{Void, AbstractCache}}
    # check span is sane
    @_checkspan(span, z)

    # init and final times
    t, T = span

    # store initial state in monitors
    _ismonitor(M) && push!(mon, t, _state(z))

    # run until condition
    while t != T
        # obtain time step from hook and what the next time is
        t_next, Δt = _next_Δt(t, T, hook(system.g, system.A, z))

        # advance
        step!(method, system, t, Δt, z, cache)

        # update
        t = t_next

        # store solution into monitor
        _ismonitor(M) && push!(mon, t, _state(z))
    end

    return z
end

# ---------------------------------------------------------------------------- #
# TIME STEPPING BASED ON CACHED STAGES, ONLY FOR LINEARISED EQUATIONS
function _propagate!(method::AbstractMethod{Z, NS, :LIN, ISADJ},
                   stepping::TimeStepFromCache{NS},
                     system::System,
                       span::_Span,
                          z::Z,
                        mon::M) where {Z, NS, ISADJ,
                                       M<:Union{Void, AbstractMonitor}}
    # check span is sane
    @_checkspan(span, z)

    # store initial state in monitors
    _ismonitor(M) && push!(mon, ts[1], _state(z))

    # TODO: fix this with proper iteration support for the stage cache
    ts  = stepping.stcache.ts
    Δts = stepping.stcache.Δts
    xs  = stepping.stcache.xs

    # integrate forward or backward based on type of linear equation
    rng = ISADJ == true ? reverse(1:length(ts)) : 1:length(ts)
    for i in rng
        step!(method, system, t[i], Δt[i], z, stages[i])
        _ismonitor(M) && push!(mon, t[i]+Δt[i], _state(z))
    end

    return z
end

# return next time and time step that allows hitting the end point exactly
function _next_Δt(t, T, Δt::S) where {S<:Real}
    @assert Δt > 0 "negative time step encountered"
    t_next = ifelse(t ≤ T, min(t+Δt, T), max(t-Δt, T))
    return t_next, S(t_next - t)
end