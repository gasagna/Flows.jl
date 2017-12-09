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
@inline (I::Integrator)(x, span::NTuple{2, Real}, mon::Vararg{Monitor}) =
    _propagate!(I.scheme, I.system, span, I.Δt, x, mon...)

# Integrator augmented with a quadrature function are callable with an additional argument.
@inline (I::Integrator)(x, q, span::NTuple{2, Real}, mon::Vararg{Monitor}) =
    _propagate!(I.scheme, I.system, span, I.Δt, aug_state(x, q), mon...)

# Main propagation function
@generated function _propagate!(scheme::IMEXRKScheme{S},
                                system::System,
                                span::NTuple{2, Real},
                                Δt::Float64, # enforced to be positive
                                z::S,
                                mon::Vararg{Monitor, N}) where {S, N}
    quote
        $(Expr(:meta, :inline))

        # define integration times
        ts = LossLessRange(span[1], span[2], span[1] < span[2] ? Δt : -Δt)

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