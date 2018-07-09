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

# Main entry points. Integrators are callable objects ...
@inline (I::Integrator)(x, span::NTuple{2, Real}, mon::Vararg{AbstractMonitor}) =
    _propagate!(I.scheme, I.system, Float64.(span), I.stepping, x, mon...)

# Integrator augmented with a quadrature function are callable with an additional argument.
@inline (I::Integrator)(x, q, span::NTuple{2, Real}, mon::Vararg{AbstractMonitor}) =
    _propagate!(I.scheme, I.system, Float64.(span), I.stepping, aug_state(x, q), mon...)

# Propagation function for constant time stepping
@generated function _propagate!(scheme::Scheme{S},
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

        while t < T
            # obtain time step from hook
            Δt = next_Δt(t, T, hook(system.g, system.A, z))

            # advance
            step!(scheme, system, t, Δt, z)

            # update
            t += Δt

            # store solution into monitor
            Base.Cartesian.@nexprs $N i->push!(mon[i], t, _state(z))
        end
        z
    end
end

# allow hitting the end point exactly
next_Δt(t, T, Δt::S) where {S} = S(min(t + Δt, T) - t)
