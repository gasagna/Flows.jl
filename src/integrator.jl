export integrator

struct Integrator{S<:AugmentedSystem, Sc}
         sys::S       # 
      scheme::Sc      # the scheme, with storage and RK implementation
          Δt::Float64 # the time step. Will be replaced by integration options
    function Integrator(sys::S,
                               scheme::Sc,
                               Δt::Real) where {S, Sc}
        Δt > 0 || throw(ArgumentError("Δt must be positive, got $Δt"))
        new{S, Sc}(sys, scheme, Δt)
    end
end

# Outer constructor
integrator(g, A, q, scheme::IMEXRKScheme, Δt::Real) = Integrator(aug_system(g, A,             q), scheme, Δt)
integrator(g, A,    scheme::IMEXRKScheme, Δt::Real) = Integrator(aug_system(g, A,       nothing), scheme, Δt)
integrator(g,       scheme::IMEXRKScheme, Δt::Real) = Integrator(aug_system(g, nothing, nothing), scheme, Δt)

# Main entry points. Integrators are callable objects ...
(I::Integrator)(x, span::NTuple{2, Real}, mon::Union{Void, Monitor}=nothing) =
    _propagate!(I.scheme, I.sys, span, I.Δt, x, mon)

# Integrator augmented with a quadrature function are callable with an additional argument.
(I::Integrator{<:AugmentedSystem})(x, q, span::NTuple{2, Real}, mon::Union{Void, Monitor}=nothing) =
    _propagate!(I.scheme, I.sys, span, I.Δt, aug_state(x, q), mon)

# Main propagation function
@inline function _propagate!(scheme::IMEXRKScheme{S}, sys, span::NTuple{2, Real}, Δt::Real, z::S, mon::Union{Void, Monitor}) where {S}
    # initial integration time
    t = Float64(span[1])

    # might wish to store initial state
    mon isa Monitor && push!(mon, t, _state(z))

    # start integration
    while integrate(t, span, Δt)
        # compute next time step
        Δt⁺ = next_Δt(t, span, Δt)

        # advance
        step!(scheme, sys, t, Δt⁺, z)

        # update time
        t += Δt⁺

        # store solution into monitor
        mon isa Monitor && push!(mon, t, _state(z))
    end
    z
end

# Evaluate condition to continue integration. This depends on the
# direction of time, which can be negative for the adjoint problem
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