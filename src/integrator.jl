using ArgCheck

export forwmap!

"""
    forwmap!(f, T, Δt, scheme)

Returns a function for the `T`-time forward map associated to the dynamical system
`f`. Integration is performed using the IMEXRK scheme defined by `scheme` using
fixed time step `Δt`. The signature of the returned function `ret` is `ret(x)`, 
which operates in place, overwriting its argument. The input argument `x` should
be of a type with the storage defined in `scheme`.
"""
function forwmap!(f::AbstractIMEXSystem, T::Real, Δt::Real, scheme::IMEXRKScheme)
    @argcheck T  > 0
    @argcheck ΔT > 0
    # the returned function will work in place and will propagate forward
    # the system by `T`
    function wrapped(x::AbstractVector)
        # set time to zero. Change this for non-autonomous systems.
        t = zero(Δt) 
        # loop
        while true
            # obtain next time step
            Δt⁺ = next_Δt(t, T, Δt)

            # step forward
            step!(scheme, nonstiff(f), stiff(f), t, Δt⁺, x)
            t += Δt⁺

            # stop when we reach final time.
            if t ≥ T 
                return x
            end
        end
    end
end

"""
    Return time step for current RK step. Becomes smaller than `Δt` in 
    case we need to hit `T` exactly.
""" 
function next_Δt(t, T, Δt)
    min(t + Δt, T) - t
end