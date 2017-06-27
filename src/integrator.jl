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
function forwmap!(g, A, Δt::Real, scheme::IMEXRKScheme)
    Δt > 0 || throw(ArgumentError("Δt must be greater than 0, got $Δt"))
    # the returned function will work in place
    function wrapped(x, T::Real, monitors::Monitor...)
        # time must be positive
        T  > 0 || throw(ArgumentError("T must be greater than 0, got $T"))

        # set time to zero. Change this for non-autonomous systems.
        t = zero(Δt) 
        
        # loop
        while true
            # obtain next time step
            Δt⁺ = next_Δt(t, T, Δt)

            # step forward
            step!(scheme, g, A, t, Δt⁺, x)
            t += Δt⁺

            # push a new sample to all monitors
            push!(monitors, x, t)

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
function next_Δt(t, T, Δt::S)::S where S
    min(t + Δt, T) - t
end