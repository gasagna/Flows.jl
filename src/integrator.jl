export integrator, fwdmapgen, _propagate!

# """
#     forwmap!(g, A, T, Δt, scheme)

# Returns a function for the `T`-time forward map associated to the dynamical system
# defined by `g` and `A`. These two define the non-stiff and stiff part of the 
# equations, and obey the interface
    
#     g(::Real, x::T, ẋ::T)
#     A_mul_B!(out::T, A, x::T)
#     ImcA!(out::T, A, x::T)

# where `T` is any custom type. The code here is agnostic to this type, as long as
# there exists a method for `similar(::Type{T})`, so that the temporaries needed 
# can be generated internally without user intervention. 

# Integration is performed using the IMEXRK scheme defined by `scheme` using
# fixed time step `Δt`. The signature of the returned function `ret` is `ret(x)`, 
# which operates in place, overwriting its argument. The input argument `x` should
# be of a type with the storage defined in `scheme`.
# """

struct Integrator{G, At, Sc}
           g::G       # the non-stiff part. Potentially augmented with a quadrature function
           A::At      # the linear stiff part.
      scheme::Sc      # the scheme, with storage and RK implementation
          Δt::Float64 # the time step. Will be replaced by integration options
    function Integrator{G, At, Sc}(g::G, 
                                   A::At, 
                                   scheme::Sc, 
                                   Δt::Real) where {G, At, Sc}
        Δt == 0 && throw(ArgumentError("Δt must be different from 0, got $Δt"))
        new(g, A, scheme, Δt)
    end
end

# Outer constructor: if a quadrature function is provided, augment system and call outer constructor
integrator(g, A,    scheme::IMEXRKScheme, Δt::Real) = Integrator{typeof.((g, A, scheme))...}(g, A, scheme, Δt)
integrator(g, A, q, scheme::IMEXRKScheme, Δt::Real) = integrator(aug_system(g, q), A, scheme, Δt)

# Main entry points. Integrators are callable objects ...
(I::Integrator)(x, T::Real)                                     = _propagate!(I.scheme, I.g, I.A, T, I.Δt, x, Work(nothing, nothing))
(I::Integrator)(x, T::Real, mon::Monitor)                       = _propagate!(I.scheme, I.g, I.A, T, I.Δt, x, Work(nothing,     mon))
(I::Integrator)(x, T::Real, sol::AbstractStorage)               = _propagate!(I.scheme, I.g, I.A, T, I.Δt, x, Work(sol,     nothing))
(I::Integrator)(x, T::Real, mon::Monitor, sol::AbstractStorage) = _propagate!(I.scheme, I.g, I.A, T, I.Δt, x, Work(mon,         sol))

# Integrator augmented with a quadrature function are callable with an additional argument.
(I::Integrator{<:AugmentedSystem})(x, q, T::Real)                                     = _propagate!(I.scheme, I.g, I.A, T, I.Δt, aug_state(x, q), Work(nothing, nothing))
(I::Integrator{<:AugmentedSystem})(x, q, T::Real, mon::Monitor)                       = _propagate!(I.scheme, I.g, I.A, T, I.Δt, aug_state(x, q), Work(nothing,     mon))
(I::Integrator{<:AugmentedSystem})(x, q, T::Real, sol::AbstractStorage)               = _propagate!(I.scheme, I.g, I.A, T, I.Δt, aug_state(x, q), Work(sol,     nothing))
(I::Integrator{<:AugmentedSystem})(x, q, T::Real, mon::Monitor, sol::AbstractStorage) = _propagate!(I.scheme, I.g, I.A, T, I.Δt, aug_state(x, q), Work(mon,         sol))

# This is used to dispatch appropriate methods
struct Work{SolType, MonType}
    sol::SolType
    mon::MonType
end

hasmonitor(w::Work)                    = true
hasmonitor(w::Work{T, Void}) where {T} = false
hasstorage(w::Work)                    = true
hasstorage(w::Work{Void, T}) where {T} = false

# Main propagation function
@inline function _propagate!(scheme::IMEXRKScheme{S}, g, A, T::Real, Δt::Real, z::S, work::Work) where {S}
    # disallow crazy stuff
    T  > 0 || throw(ArgumentError("T must be greater than 0, got $T"))

    # initial integration time depends on integration mode
    t = Δt > 0 ? zero(T) : T 

    # might wish to store initial state in the storage and monitor
    hasstorage(work) && iswmode(work.sol) && push!(work.sol, _state(z), t)
    hasmonitor(work) && push!(work.mon, t, _state_quad(z))

    # start integration
    while integrate(0, t, T, Δt)
        # compute next time step
        Δt⁺ = next_Δt(0, t, T, Δt)

        # advance (might need storage, for reading or storing)
        hasstorage(work) && isrmode(work.sol) && (step!(scheme, g, A, t, Δt⁺, z, storage))
        hasstorage(work) && iswmode(work.sol) && (step!(scheme, g, A, t, Δt⁺, z); push!(work.sol, _state(z), t+Δt⁺))
        hasstorage(work) || step!(scheme, g, A, t, Δt⁺, z)

        # update monitor if needed
        hasmonitor(work) && push!(ms, t, _state_quad(z))

        # update time
        t += Δt⁺
    end
    z
end

# Evaluate condition to continue integration. This depends on the
# direction of time, which can be negative for the adjoint problem
function integrate(tmin, t, tmax, Δt)
    Δt > 0 && return t < tmax # forward
    Δt < 0 && return t > tmin # backward
end

# Return time step for current RK step. Becomes smaller than `Δt` in 
# case we need to hit the stopping `T` exactly.
function next_Δt(tmin, t, tmax, Δt::S)::S where S
    Δt > 0 && return min(t + Δt, tmax) - t # forward
    Δt < 0 && return max(t + Δt, tmin) - t # backward
end