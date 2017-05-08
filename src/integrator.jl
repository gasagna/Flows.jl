using ArgCheck

export forwmap!

"""
    forwmap!(g, A, T, Δt, scheme)

Returns a function for the `T`-time forward map associated to the dynamical system
defined by `g` and `A`. These two define the non-stiff and stiff part of the 
equations, and obey the interface
    
    g(::Real, x::T, ẋ::T)
    A_mul_B!(out::T, A, x::T)
    ImcA!(out::T, A, x::T)

where `T` is any custom type. The code here is agnostic to this type, as long as
there exists a method for `similar(::Type{T})`, so that the temporaries needed 
can be generated internally without user intervention. 

Integration is performed using the IMEXRK scheme defined by `scheme` using
fixed time step `Δt`. The signature of the returned function `ret` is `ret(x)`, 
which operates in place, overwriting its argument. The input argument `x` should
be of a type with the storage defined in `scheme`.
"""
function forwmap!(g, A, T::Real, Δt::Real, scheme::IMEXRKScheme)
    @argcheck T  > 0
    @argcheck Δt > 0
    # the returned function will work in place and will propagate forward
    # the system by `T`
    function wrapped(x)
        # set time to zero. Change this for non-autonomous systems.
        t = zero(Δt) 
        # loop
        while true
            # obtain next time step
            Δt⁺ = next_Δt(t, T, Δt)

            # step forward
            step!(scheme, g, A, t, Δt⁺, x)
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