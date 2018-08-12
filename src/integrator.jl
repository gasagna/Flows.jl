export flow, InvalidSpanError

# ---------------------------------------------------------------------------- #
# The Flow type!
struct Flow{TS<:AbstractTimeStepping, M<:AbstractMethod, S<:System}
    tstep::TS # the method used for time stepping
     meth::M  # the method, with storage, implementation and time stepping
      sys::S  # the system to be integrated
    Flow(ts::TS, m::M, sys::S) where {TS, M, S} = new{TS, M, S}(ts, m, sys)
end

# main entry point specifying explicit, implicit and quadrature parts
flow(g, A, q, m::AbstractMethod, ts::AbstractTimeStepping) =
     Flow(ts, m, System(g, A, q))

flow(g, A,    m::AbstractMethod, ts::AbstractTimeStepping) =
     Flow(ts, m, System(g, A, nothing))

flow(g,       m::AbstractMethod, ts::AbstractTimeStepping) =
     Flow(ts, m, System(g, nothing, nothing))

# ---------------------------------------------------------------------------- #
# FLOWS ARE CALLABLE OBJECTS: THIS IS THE MAIN INTERFACE

# normal stepping
(I::Flow)(x, span::NTuple{2, Real}) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span), x, nothing, nothing)

# fill a monitor
(I::Flow)(x, span::NTuple{2, Real}, m::AbstractMonitor) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span), x, nothing, m)

# fill a cache and optionally a monitor
(I::Flow)(x, span::NTuple{2, Real}, 
          c::AbstractStageCache, m::Union{Void, <:AbstractMonitor}=nothing) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span), x, c, m)

# stepping based on cache only 
(I::Flow{TimeStepFromCache})(x, 
                             c::AbstractStageCache, 
                             m::Union{Void, <:AbstractMonitor}=nothing) =
    _propagate!(I.meth, I.sys, x, c, m)

# ---------------------------------------------------------------------------- #
# SAME WITH QUADRATURE
# normal stepping
(I::Flow)(x, q, span::NTuple{2, Real}) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span), 
                coupled(x, q), nothing, nothing)

# fill a monitor
(I::Flow)(x, q, span::NTuple{2, Real}, m::AbstractMonitor) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span), 
                coupled(x, q), nothing, m)

# fill a cache and optionally a monitor
(I::Flow)(x, q, span::NTuple{2, Real}, 
          c::AbstractStageCache, m::Union{Void, <:AbstractMonitor}=nothing) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span), coupled(x, q), c, m)

# stepping based on cache only, calculating a quadrature
(I::Flow{TimeStepFromCache})(x, q,
                             c::AbstractStageCache, 
                             m::Union{Void, <:AbstractMonitor}=nothing) =
    _propagate!(I.meth, I.sys, coupled(x, q), c, m)


# ---------------------------------------------------------------------------- #
# PROPAGATION FUNCTIONS & UTILS

struct InvalidSpanError{S} <: Exception 
    span::S
end

Base.showerror(io::IO, e::InvalidSpanError) =
    print(io, "Invalid time span ", e.span, ". Time must be increasing.\n")

# check span and/or return
macro _checkspan(span, z)
    quote
       $(esc(span))[1] ==$(esc(span))[2] && return $z
       $(esc(span))[1]  >$(esc(span))[2] && throw(InvalidSpanError($(esc(span))))
    end
end

# ---------------------------------------------------------------------------- #
# CONSTANT FORWARD TIME STEPPING, ONLY FOR STATE EQUATIONS
function _propagate!(method::AbstractMethod{Z, NS, :NORMAL},
                   stepping::TimeStepConstant,
                     system::System,
                       span::NTuple{2, Real},
                          z::Z,
                      cache::C,
                        mon::M) where {Z, NS,
                                       M<:Union{Void, AbstractMonitor},
                                       C<:Union{Void, AbstractStageCache}}
    # check span is sane
    @_checkspan(span, z)

    # define integration times
    ts = LossLessRange(span[1], span[2], stepping.Δt)

    # push initial state monitor
    _ismonitor(M) && push!(mon, ts[1], first(z))

    # start integration
    for j = 2:length(ts)
        step!(method, system, ts[j-1], ts[j]-ts[j-1], z, cache)
        _ismonitor(M) && push!(mon, ts[j], first(z))
    end

    return z
end

# ---------------------------------------------------------------------------- #
# PROPAGATION BASED ON SYSTEMS HOOK: ONLY FOR STATE EQUATIONS
function _propagate!(method::AbstractMethod{Z, NS, :NORMAL},
                       hook::AbstractTimeStepFromHook,
                     system::System,
                       span::NTuple{2, Real},
                          z::Z,
                      cache::C,
                        mon::M) where {Z, NS,
                                       M<:Union{Void, AbstractMonitor},
                                       C<:Union{Void, AbstractStageCache}}
    # check span is sane
    @_checkspan(span, z)

    # init and final times
    t, T = span

    # store initial state in monitors
    _ismonitor(M) && push!(mon, t, first(z))

    # run until condition
    while t != T
        # obtain time step from hook and what the next time is
        t_next, Δt = _next_Δt(t, T, hook(system.g, system.A, z))

        # advance
        step!(method, system, t, Δt, z, cache)

        # update
        t = t_next

        # store solution into monitor
        _ismonitor(M) && push!(mon, t, first(z))
    end

    return z
end

# return next time and time step that allows hitting the end point exactly
function _next_Δt(t, T, Δt::S) where {S<:Real}
    @assert Δt > 0 "negative time step encountered"
    t_next = ifelse(t ≤ T, min(t+Δt, T), max(t-Δt, T))
    return t_next, S(t_next - t)
end

# ---------------------------------------------------------------------------- #
# TIME STEPPING BASED ON CACHED STAGES, ONLY FOR LINEARISED EQUATIONS
function _propagate!(method::AbstractMethod{Z, NS, :LIN, ISADJ},
                     system::System,
                          z::Z,
                      cache::AbstractStageCache{NS},
                        mon::M) where {Z, NS, ISADJ,
                                       M<:Union{Void, AbstractMonitor}}
    # store initial state in monitors
    _ismonitor(M) && push!(mon, ts[1], first(z))

    # TODO: fix this with proper iteration support for the stage cache
    ts  = cache.ts
    Δts = cache.Δts
    xs  = cache.xs

    # integrate forward or backward based on type of linear equation
    rng = ISADJ == true ? reverse(1:length(ts)) : 1:length(ts)
    for i in rng
        step!(method, system, ts[i], Δts[i], z, xs[i])
        _ismonitor(M) && push!(mon, ts[i]+Δts[i], first(z))
    end

    return z
end