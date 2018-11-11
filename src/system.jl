import LinearAlgebra: mul!

# ~~~ TYPE FOR INTEGRATING N COUPLED PROBLEMS ~~~
struct System{N, GT, AT}
    g::GT # explicit term
    A::AT # linear implicit term
end

# constructors
System(g::GT, A::AT) where {GT, AT} = System{1, GT, AT}(g, A)
System(g::GT, A::AT) where {N, GT<:Coupled{N}, AT<:Coupled{N}} = 
    System{N, GT, AT}(g, A)

# ~ Explicit part ~
(sys::System{1})(t::Real, z, dzdt) = sys.g(t, z, dzdt)

(sys::System{2})(t::Real, z::Coupled{2}, dzdt::Coupled{2}) =
    (sys.g[1](t, z[1], dzdt[1]);
     sys.g[2](t, z[1], dzdt[1], z[2], dzdt[2]); dzdt)

(sys::System{3})(t::Real, z::Coupled{3}, dzdt::Coupled{3}) =
     (sys.g[1](t, z[1], dzdt[1]);
      sys.g[2](t, z[1], dzdt[1], z[2], dzdt[2]);
      sys.g[3](t, z[1], dzdt[1], z[2], dzdt[2], z[3], dzdt[3]); dzdt)

# add generic method for larger N if needed

# this is for the adjoint schemes
(sys::System{1})(t::Real, u, z, dzdt) = sys.g(t, u, z, dzdt)

# ~ Implicit part I ~
mul!(out, sys::System{1, GT, AT}, z) where {GT, AT} =
    ((AT isa Nothing ? (out .= 0) : mul!(out, sys.A, z)); out)

@generated function mul!(out::Coupled{N},
                         sys::System{N, Coupled{N, GT}, Coupled{N, AT}},
                           z::Coupled{N}) where {N, GT, AT}
    return quote
        Base.Cartesian.@nexprs $N i->($(AT.parameters)[i] == Nothing ?
                (out[i] .= 0) : mul!(out[i], sys.A[i], z[i]))
        return out
    end
end

# ~ Implicit part II ~
# Since we treat the quadrature fully explicitly, the solution of
# (I-cA)z = y for the quadrature part is simply z = y, because the
# component of A associated to this part is zero and the state
# and quadrature parts are decoupled.
ImcA!(sys::System{1, GT, AT}, c::Real, y, z) where {GT, AT} =
    ((AT isa Nothing ? (z .= y) : ImcA!(sys.A, c, y, z)); z)

@generated function ImcA!(sys::System{N, Coupled{N, GT}, Coupled{N, AT}},
                            c::Real,
                            y::Coupled{N},
                            z::Coupled{N}) where {N, GT, AT}
    return quote 
        Base.Cartesian.@nexprs $N i->($(AT.parameters)[i] == Nothing ?
                (z[i] .= y[i]) : ImcA!(sys.A[i], c, y[i], z[i]))
        return z
    end
end