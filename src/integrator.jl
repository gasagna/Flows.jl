using ArgCheck

export forwmap!

function forwmap!(f::AbstractIMEXSystem, T::Real, Δt::Real, impl::IMEXRKScheme)
    @argcheck T > 0
    # the returned function will work in place and will propagate forward
    # the system by `T`
    function wrapped(x::AbstractVector)
        # set time to zero. Change this for non-autonomous systems.
        t = zero(T)
        # loop
        while true
            # obtain next time step
            Δt⁺ = next_Δt(t, T, Δt)

            # step forward
            step!(impl, nonstiff(f), stiff(f), t, Δt⁺, x)
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