# Type for ODE/PDE problems
struct System{G, A, Q}
    g::G # explicit term
    A::A # linear implicit term: can be Void if not provided
    q::Q # quadrature function: can be Void if not provided
end
System(g, A, q) = System(g, A, q)

function (aug::System{G, A, Q})(t::Real, z::T, dzdt::T) where {G, A, Q, T} 
    aug.g(t, _state(z), _state(dzdt))
    Q <: Void || aug.q(t, _state(z), _quad(dzdt))
    return dzdt
end

function A_mul_B!(out::T, aug::System{G, A}, z::T) where {G, A, T} 
    if A <: Void
        out .*= 0
    else
        A_mul_B!(_state(out), aug.A, _state(z))
    end
    return out
end

# Since we treat the quadrature fully explicitly, the solution of
# (I-cA)z = y for the quadrature part is simply z = y, because the
# component of A associated to this part is zero and the state 
# and quadrature parts are decoupled.
function ImcA!(aug::System{G, A, Q}, c::Real, y::T, z::T) where {G, A, Q, T} 
    if A <: Void
        z .= y
    else
        ImcA!(aug.A, c, _state(y), _state(z))
        T<:AugmentedState && _quad(z) .= _quad(y)
    end
    return z
end