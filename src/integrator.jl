export integrator

struct Integrator{S<:System, Sc}
      system::S       # the system to be integrated
      scheme::Sc      # the scheme, with storage and RK implementation
          Δt::Float64 # the time step. Will be replaced by integration options
    function Integrator(system::S,
                        scheme::Sc,
                        Δt::Real) where {S, Sc}
        Δt > 0 || throw(ArgumentError("Δt must be positive, got $Δt"))
        new{S, Sc}(system, scheme, Δt)
    end
end

# Outer constructor
integrator(g, A, q, scheme::IMEXRKScheme, Δt::Real) = Integrator(System(g, A,             q), scheme, Δt)
integrator(g, A,    scheme::IMEXRKScheme, Δt::Real) = Integrator(System(g, A,       nothing), scheme, Δt)
integrator(g,       scheme::IMEXRKScheme, Δt::Real) = Integrator(System(g, nothing, nothing), scheme, Δt)

# Main entry points. Integrators are callable objects ...
(I::Integrator)(x, span::NTuple{2, Real}, mon::Vararg{Monitor}) =
    _propagate!(I.scheme, I.system, span, I.Δt, x, mon...)

# Integrator augmented with a quadrature function are callable with an additional argument.
(I::Integrator)(x, q, span::NTuple{2, Real}, mon::Vararg{Monitor}) =
    _propagate!(I.scheme, I.system, span, I.Δt, aug_state(x, q), mon...)

# Main propagation function
@generated function _propagate!(scheme::IMEXRKScheme{S}, system::System, span::NTuple{2, Real}, Δt::Real, z::S, mon::Vararg{Monitor, N}) where {S, N}
    quote
        $(Expr(:meta, :inline))
        
        # initial integration time
        t = Float64(span[1])

        # might wish to store initial state in monitors
        $N > 0 && (Base.Cartesian.@nexprs $N i->push!(mon[i], t, _state(z)))

        # start integration
        while integrate(t, span, Δt)
            # compute next time step
            Δt⁺ = next_Δt(t, span, Δt)

            # advance
            step!(scheme, system, t, Δt⁺, z)

            # update time
            t += Δt⁺

            # store solution into monitor
            $N > 0 && (Base.Cartesian.@nexprs $N i->push!(mon[i], t, _state(z)))
        end
        z
    end
end

# Evaluate condition to continue integration. This depends on the
# direction of time, which can be negative for the adjoint system
function integrate(t, span, Δt)
    tstart, tstop = span 
    tstart ≤ tstop && return t < tstop # forward
    tstart > tstop && return t > tstop # backward
end

# Return time step for current step. Becomes smaller than
# `Δt` in case we need to hit the stopping time exactly
function next_Δt(t, span, Δt::S)::S where S
    tstart, tstop = span 
    tstart ≤ tstop && return min(t + Δt, tstop) - t # forward
    tstart > tstop && return max(t - Δt, tstop) - t # backward
end