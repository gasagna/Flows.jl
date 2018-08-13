# Type for ODE/PDE problems
struct System{G, A, Q}
    g::G # explicit term
    A::A # linear implicit term: can be Void if not provided
    q::Q # quadrature function: can be Void if not provided
end

# Explicit part
(sys::System{G, A, Q})(t::Real, z, dzdt) where {G, A, Q} =
    (sys.g(t, first(z), first(dzdt)); 
     sys.q(t, first(z),  last(dzdt)); dzdt)

(sys::System{G, A, Void})(t::Real, z, dzdt) where {G, A} = 
    (sys.g(t, z, dzdt); dzdt)

# Explicit part for linearised problems
(sys::System{G, A, Q})(t::Real, u, z, dzdt) where {G, A, Q} =
    (sys.g(t, u, first(z), first(dzdt));
     sys.q(t,    first(z),  last(dzdt)); dzdt)

(sys::System{G, A, Void})(t::Real, u, z, dzdt) where {G, A} = 
    (sys.g(t, u, z, dzdt); dzdt)

# Explicit part for coupled integration. We execute the FIRST function first,
# then use the output for evaluating the second function too. This assumes
# essentially that the coupled system is lower triangular. If passing a
# quadrature, we expect if will work with Coupled objects. The way this is 
# coded up does not assume that the implicit part is coupled too. Note that
# these two will not get called if only `z` and `dzdt` are `Coupled`. This means
# that for quadrature integration the previous methods will apply. This is
# a bit auto-magical, but works well.
# TODO: we might want to pass dzdt to the quadrature too, if needed.
function (sys::System{<:Coupled, A, Q})(t::Real,
                                       z::Coupled{C},
                                    dzdt::Coupled{C}) where {A, Q, C<:Coupled}
    # unpack things
     x, y,   q = first(first(z)),    last(first(z)),    last(z)
    dx, dy, dq = first(first(dzdt)), last(first(dzdt)), last(dzdt)

    # call nonlinear equations first
    first(sys.g)(t, x, dx)

    # call linearised equations second
    last(sys.g)(t, x, dx, y, dy)

    # call quadrature
    sys.q(t, first(z), dq)

    # return
    return dzdt
end

(sys::System{<:Coupled, A, Void})(t::Real, z::Coupled, dzdt::Coupled) where {A} =
    (first(sys.g)(t, first(z), first(dzdt));
      last(sys.g)(t, first(z), first(dzdt), last(z), last(dzdt)); dzdt)


# Implicit part. We also define methods for At_mul_B!, for the adjoint
# code, hence the user must also provide At_mul_B! for his linear implicit
# type. Note we set the quadrature part to zero, since we are treating it
# fully explicitly.
for fun in [:A_mul_B!, :At_mul_B!]
    @eval begin
        Base.$fun(out, sys::System{G, A, Q}, z) where {G, A, Q} =
            ($fun(first(out), sys.A, first(z)); last(out) .= 0; out)

        Base.$fun(out, sys::System{G, A, Void}, z) where {G, A} =
            $fun(out, sys.A, z)

        # when A is a `Coupled` object. Obviously, we need `y` and `z` to be
        # `Coupled` as well - this is with quadrature
        Base.$fun(out::Coupled{Coupled},
                  sys::System{G, <:Coupled, Q},
                    z::Coupled{Coupled}) where {G, Q} =
            ($fun(first(first(out)), first(sys.A), first(first(z)));
             $fun(first( last(out)),  last(sys.A), first( last(z)));
             last(out) .= 0; out)

        # and for no quadrature
        Base.$fun(out::Coupled, 
                  sys::System{G, <:Coupled, Void}, 
                    z::Coupled) where {G} =
            ($fun(first(out), first(sys.A), first(z));
             $fun( last(out),  last(sys.A),  last(z)); out)

        # this is for fully explicit problems. We do not add methods for when
        # G, out and z are `Coupled` objects, since `out::Coupled` should
        # support broadcasting operations
        Base.$fun(out, sys::System{G, Void, Q}, z) where {G, Q} =
            (out .= 0; out)

        Base.$fun(out, sys::System{G, Void, Void}, z) where {G} = 
            (out .= 0; out)
    end
end

# Since we treat the quadrature fully explicitly, the solution of
# (I-cA)z = y for the quadrature part is simply z = y, because the
# component of A associated to this part is zero and the state
# and quadrature parts are decoupled. Same for no linear term.
for fun in [:ImcA!, :ImcAt!]
    @eval begin
        $fun(sys::System{G, A, Q}, c::Real, y, z) where {G, A, Q} =
            ($fun(sys.A, c, first(y), first(z)); last(z) .= last(y); z)

        $fun(sys::System{G, A, Void}, c::Real, y, z) where {G, A} =
            ($fun(sys.A, c, y, z); z)

        # when A is a `Coupled` object. Obviously, we need `y` and `z` to be
        # `Coupled` as well - this is with quadrature.
        $fun(sys::System{G, <:Coupled, Q}, 
               c::Real, 
               y::Coupled{Coupled}, 
               z::Coupled{Coupled}) where {G, Q} =
            ($fun(first(sys.A), c, first(first(y)), first(first(z)));
             $fun( last(sys.A), c, first( last(y)), first( last(z)));
                last(z) .= last(y); z)

        # and with no quadrature.
        $fun(sys::System{G, <:Coupled, Void},
               c::Real,
               y::Coupled,
               z::Coupled) where {G} =
            ($fun(first(sys.A), c, first(y), first(z));
             $fun( last(sys.A), c,  last(y),  last(z)); z)

        # this is for fully explicit systems. We do not add methods for when
        # G, y and z are `Coupled` objects, since `z` and `y` as `Coupled` objects 
        # should support broadcasting operations
        $fun(sys::System{G, Void, Q}, c::Real, y, z) where {G, Q} =
            (z .= y; z)

        $fun(sys::System{G, Void, Void}, c::Real, y, z) where {G} =
            (z .= y; z)
    end
end