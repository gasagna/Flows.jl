# Type for ODE/PDE problems
struct System{G, A, Q}
    g::G # explicit term
    A::A # linear implicit term: can be Void if not provided
    q::Q # quadrature function: can be Void if not provided
end

# Explicit part
(sys::System{G, A, Q})(t::Real, z, dzdt) where {G, A, Q} =
    (sys.g(t, _state(z), _state(dzdt)); sys.q(t, _state(z), _quad(dzdt)); dzdt)

(sys::System{G, A, Void})(t::Real, z, dzdt) where {G, A} = 
    (sys.g(t, z, dzdt); dzdt)

# Implicit part
Base.A_mul_B!(out, sys::System{G, A, Q}, z) where {G, A, Q} =
    (A_mul_B!(_state(out), sys.A, _state(z)); _quad(out) .= 0; out)

Base.A_mul_B!(out, sys::System{G, A, Void}, z) where {G, A} =
    A_mul_B!(out, sys.A, z)

Base.A_mul_B!(out, sys::System{G, Void, Q}, z) where {G, Q} =
    (out .= 0; out)

Base.A_mul_B!(out, sys::System{G, Void, Void}, z) where {G} = 
    (out .= 0; out)

# Since we treat the quadrature fully explicitly, the solution of
# (I-cA)z = y for the quadrature part is simply z = y, because the
# component of A associated to this part is zero and the state
# and quadrature parts are decoupled. Same for no linear term.
ImcA!(sys::System{G, A, Q}, c::Real, y, z) where {G, A, Q} =
    (ImcA!(sys.A, c, _state(y), _state(z)); _quad(z) .= _quad(y); z)

ImcA!(sys::System{G, A, Void}, c::Real, y, z) where {G, A} =
    (ImcA!(sys.A, c, y, z); z)

ImcA!(sys::System{G, Void, Q}, c::Real, y, z) where {G, Q} =
    (z .= y; z)

ImcA!(sys::System{G, Void, Void}, c::Real, y, z) where {G} =
    (z .= y; z)