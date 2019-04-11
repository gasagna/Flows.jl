export flow, InvalidSpanError

# ---------------------------------------------------------------------------- #
# The Flow type!
mutable struct Flow{TS<:AbstractTimeStepping, M<:AbstractMethod, S<:System}
    tstep::TS # the method used for time stepping
     meth::M  # the method, with storage, implementation and time stepping
      sys::S  # the system to be integrated
    Flow(ts::TS, m::M, sys::S) where {TS, M, S} = new{TS, M, S}(ts, m, sys)
end

# single system
flow(g, m::AbstractMethod, ts::AbstractTimeStepping) =
    flow(g, nothing, m, ts)

flow(g, A, m::AbstractMethod, ts::AbstractTimeStepping) =
     Flow(ts, m, System(g, A))

 # coupled systems
flow(g::Coupled{N}, m::AbstractMethod, ts::AbstractTimeStepping) where {N} =
     flow(g, default_dep(N), m, ts)

flow(g::Coupled{N}, spec::CallDependency{N}, m::AbstractMethod, 
                                            ts::AbstractTimeStepping) where {N} =
     flow(g, couple(ntuple(i->nothing, N)...), spec, m, ts)

flow(g::Coupled{N}, A::Coupled{N}, m::AbstractMethod,
                                  ts::AbstractTimeStepping) where {N} =
    flow(g, A, default_dep(N), m, ts)

flow(g::Coupled{N}, A::Coupled{N}, spec::CallDependency{N},
                                  m::AbstractMethod,
                                  ts::AbstractTimeStepping) where {N} =
     Flow(ts, m, System(g, A, spec))

# ---------------------------------------------------------------------------- #
# FLOWS ARE CALLABLE OBJECTS: THIS IS THE MAIN INTERFACE

# normal stepping
(I::Flow)(x, span::NTuple{2, Real}) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span), x, nothing, nothing, nothing)

# fill a monitor
(I::Flow)(x, span::NTuple{2, Real}, m::AbstractMonitor) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span), x, nothing, nothing, m)

# fill a cache and optionally a monitor
(I::Flow)(x, span::NTuple{2, Real}, 
          c::AbstractStageCache, m::Union{Nothing, <:AbstractMonitor}=nothing) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span), x, c, nothing, m)

# fill a store and optionally a monitor
(I::Flow)(x, span::NTuple{2, Real}, 
          s::AbstractStorage, m::Union{Nothing, <:AbstractMonitor}=nothing) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span), x, nothing, s, m)

# stepping based on cache only, with optional monitor
(I::Flow{TimeStepFromCache})(x, 
                             c::AbstractStageCache, 
                             m::Union{Nothing, <:AbstractMonitor}=nothing) =
    _propagate!(I.meth, I.sys, x, c, m)

# stepping based on store over a given span with optional monitor
(I::Flow{TimeStepFromStorage})(x,
                               s::AbstractStorage,
                               span::NTuple{2, Real},
                               m::Union{Nothing, <:AbstractMonitor}=nothing) =
    _propagate!(I.meth, I.tstep, I.sys, Float64.(span), x, s, m)

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
       $(esc(span))[1] ==$(esc(span))[2] && return $(esc(z))
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
                      store::S,
                        mon::M) where {Z, NS,
                                       S<:Union{Nothing, AbstractStorage},
                                       M<:Union{Nothing, AbstractMonitor},
                                       C<:Union{Nothing, AbstractStageCache}}
    # check span is sane
    @_checkspan(span, z)

    # define integration times
    ts = LossLessRange(span[1], span[2], stepping.Δt)
    Nsteps = length(ts)

    # push initial state to monitor and storage
    _ismonitor(M) && push!(mon,   ts[1], z)
    _isstorage(S) && push!(store, ts[1], copy(z))

    # start integration
    for j = 2:Nsteps
        step!(method, system, ts[j-1], ts[j]-ts[j-1], z, cache)
        if _ismonitor(M) 
            # skip all pushes except the last but one
            M <: StoreOneButLast && (j != Nsteps - 2 && continue) 
            push!(mon, ts[j], z)
        end
        if _isstorage(S) 
            push!(store, ts[j], copy(z))
        end
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
                      store::S,
                        mon::M) where {Z, NS,
                                       S<:Union{Nothing, AbstractStorage},
                                       M<:Union{Nothing, AbstractMonitor},
                                       C<:Union{Nothing, AbstractStageCache}}
    # check span is sane
    @_checkspan(span, z)

    # init and final times
    t, T = span

    # store initial state in monitors
    _ismonitor(M) && push!(mon,   t, z)
    _ismonitor(S) && push!(store, t, copy(z))

    # run until condition
    while t != T
        # obtain time step from hook and what the next time is
        t_next, Δt = _next_Δt(t, T, hook(system.g, system.A, z))

        # advance
        step!(method, system, t, Δt, z, cache)

        # update
        t = t_next

        # store solution into monitor
        _ismonitor(M) && push!(mon, t, z)
        _ismonitor(S) && push!(store, t, copy(z))
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
                                       M<:Union{Nothing, AbstractMonitor}}
    # TODO: fix this with proper iteration support for the stage cache
    ts  = cache.ts
    Δts = cache.Δts
    xs  = cache.xs

    # integrate forward or backward based on type of linear equation
    if ISADJ == false
        # store final state in monitors. Note cache does not contain final T.
        _ismonitor(M) && push!(mon, ts[1], z)

        for i in 1:length(ts)
            # make step
            step!(method, system, ts[i], Δts[i], z, xs[i])

            # then save current state
            _ismonitor(M) && push!(mon, ts[i]+Δts[i], z)
        end
    else 
        # store final state in monitors. Note cache does not contain final T.
        _ismonitor(M) && push!(mon, ts[end] + Δts[end], z)

        for i in reverse(1:length(ts))
            step!(method, system, ts[i], Δts[i], z, xs[i])
            _ismonitor(M) && push!(mon, ts[i], z)
        end
    end

    return z
end

# ---------------------------------------------------------------------------- #
# TIME STEPPING BASED ON STORAGE FOR CONTINUOS ADJOINT/TANGENT EQUATIONS
function _propagate!(method::AbstractMethod{Z},
                   stepping::TimeStepFromStorage,
                     system::System,
                       span::NTuple{2, Any},
                          z::Z,
                      store::AbstractStorage,
                        mon::M) where {Z, M<:Union{Nothing, AbstractMonitor}}

    # Define integration times. Note that the adjoint case, where 
    # span[1] > span[2], is handled automatically here by this "ts" object
    ts = LossLessRange(span[1], span[2], stepping.Δt)
    
    # store initial state in monitors (this could be the final adjoint state)
    _ismonitor(M) && push!(mon, ts[1], z)
    
    # march in time
    for j = 2:length(ts)
        # exec step
        step!(method, system, ts[j-1], stepping.Δt, z, store)

        # store
        _ismonitor(M) && push!(mon, ts[j], z)
    end

    return z
end