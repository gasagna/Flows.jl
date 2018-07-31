# ---------------------------------------------------------------------------- #
export CB3R2R,
       CB3R2R_TAN,
       CB3R2R_ADJ

# ---------------------------------------------------------------------------- #
# Allowed tags for 3R2R schemes
const __allowed_3R2R__ = Dict(:_2 => CB2, :_3c =>CB3c, :_3e =>CB3e)

# ---------------------------------------------------------------------------- #
# Three-register version of [2R] IMEXRK schemes from section 1.2.1 of CB 2015
struct CB3R2R{X, ISCACHE, NS, TAB} <: AbstractMethod
    store::NTuple{3, X} # these are low storage methods!
    tableau::TAB
    function CB3R2R{ISCACHE}(st::NTuple{3, X}, t::TAB) where {ISCACHE, X, TAB}
        # convert tableau to use required element type. Note that when
        # `x` is a `AugmentedState` the `eltype` of the tableau is the
        # promotion of the types of state and quadrature components.
        T = promote_type(eltype.(_state_quad(st[1]))...)
        _t = convert(IMEXTableau{T, nstages(t)}, t)
        new{X, ISCACHE, nstages(_t), typeof(_t)}(st, _t)
    end
end

# outer constructor
function CB3R2R(x::X, tag::Symbol, iscache::Bool=false) where {X}
    tag ∈ keys(__allowed_3R2R__) || error("invalid scheme tag")
    return CB3R2R{iscache}(ntuple(i->similar(x), 3), __allowed_3R2R__[tag])
end

# with quadrature part provided
CB3R2R(x::X, q::Q, tag::Symbol, iscache::Bool=false) where {X, Q} =
    CB3R2R(augment(x, q), tag, iscache)

# perform step
function step!(method::CB3R2R{X, ISCACHE, NS},
                  sys::System,
                    t::Real,
                   Δt::Real,
                    x::X,
              stcache::C) where {X, ISCACHE, NS, C<:AbstractStageCache{NS, X}}

    # hoist temporaries out
    y, z, w  = method.store
    tab = method.tableau

    # temporary vector for storing stages
    ISCACHE && (stages = sizehint!(X[], NS))

    # loop over stages
    for k = 1:NS
        if k == 1
            y .= x
        else
            y .= x .+ (tab[:aᴵ, k, k-1] .- tab[:bᴵ, k-1]).*Δt.*z .+
                      (tab[:aᴱ, k, k-1] .- tab[:bᴱ, k-1]).*Δt.*y
        end
        A_mul_B!(z, sys, y)                 # compute z = A*y then
        ImcA!(sys, tab[:aᴵ, k, k]*Δt, z, z) # get z = (I-cA)⁻¹*(A*y) in place
        w .= y .+ tab[:aᴵ, k, k].*Δt.*z     # w is the temp input, output is y
        sys(t + tab[:cᴱ, k]*Δt, w, y); ISCACHE && (push!(stages, copy(w)))
        x .= x .+ tab[:bᴵ, k].*Δt.*z .+ tab[:bᴱ, k].*Δt.*y
    end

    # cache stages if requested
    ISCACHE && push!(stcache, t, Δt, tuple(stages...))

    return nothing
end


# ---------------------------------------------------------------------------- #
# Linearisation
struct CB3R2R_TAN{X, NS, TAB} <: AbstractMethod
    store::NTuple{3, X} # these are low storage methods!
    tableau::TAB
    function CB3R2R_TAN(st::NTuple{3, X}, t::TAB) where {X, TAB}
        T = promote_type(eltype.(_state_quad(st[1]))...)
        _t = convert(IMEXTableau{T, nstages(t)}, t)
        new{X, nstages(_t), typeof(_t)}(st, _t)
    end
end

# outer constructor
function CB3R2R_TAN(x::X, tag::Symbol) where {X}
    tag ∈ keys(__allowed_3R2R__) || error("invalid scheme tag")
    return CB3R2R_TAN(ntuple(i->similar(x), 3), __allowed_3R2R__[tag])
end

# with quadrature part provided
CB3R2R_TAN(x::X, q::Q, tag::Symbol, iscache::Bool=false) where {X, Q} =
    CB3R2R_TAN(augment(x, q), tag, iscache)

# perform step
function step!(method::CB3R2R_TAN{X, NS},
                  sys::System,
                    t::Real,
                   Δt::Real,
                    x::X,
               stages::NTuple{NS, X}) where {X, NS}
    # hoist temporaries out
    y, z, w  = method.store
    tab = method.tableau

    # loop over stages
    for k = 1:NS
        # F
        if k == 1
            y .= x
        else
            y .= x .+ (tab[:aᴵ, k, k-1] .- tab[:bᴵ, k-1]).*Δt.*z .+
                      (tab[:aᴱ, k, k-1] .- tab[:bᴱ, k-1]).*Δt.*y
        end
        # E
        A_mul_B!(z, sys, y)                 # compute z = A*y then
        # D
        ImcA!(sys, tab[:aᴵ, k, k]*Δt, z, z) # get z = (I-cA)⁻¹*(A*y) in place
        # C
        w .= y .+ tab[:aᴵ, k, k].*Δt.*z     # w is the temp input, output is y
        # B
        sys(t + tab[:cᴱ, k]*Δt, stages[k], w, y)
        # A
        x .= x .+ tab[:bᴵ, k].*Δt.*z .+ tab[:bᴱ, k].*Δt.*y
    end

    return nothing
end


# ---------------------------------------------------------------------------- #
# Adjoint version
struct CB3R2R_ADJ{X, NS, TAB} <: AbstractMethod
    store::NTuple{5, X} # these are low storage methods!
    tableau::TAB
    function CB3R2R_ADJ(st::NTuple{5, X}, t::TAB) where {X, TAB}
        T = promote_type(eltype.(_state_quad(st[1]))...)
        _t = convert(IMEXTableau{T, nstages(t)}, t)
        new{X, nstages(_t), typeof(_t)}(st, _t)
    end
end

# outer constructor
function CB3R2R_ADJ(x::X, tag::Symbol) where {X}
    tag ∈ keys(__allowed_3R2R__) || error("invalid scheme tag")
    return CB3R2R_ADJ(ntuple(i->similar(x), 5), __allowed_3R2R__[tag])
end

# with quadrature part provided
CB3R2R_ADJ(x::X, q::Q, tag::Symbol, iscache::Bool=false) where {X, Q} =
    CB3R2R_ADJ(augment(x, q), tag, iscache)

# perform step
function step!(method::CB3R2R_ADJ{X, NS},
                  sys::System,
                    t::Real,
                   Δt::Real,
                    x::X,
               stages::NTuple{NS, X}) where {X, NS}
    # hoist temporaries out
    y, z, w, v, s  = method.store
    y .= 0; z .= 0; w .= 0; v .= 0; s .= 0;
    tab = method.tableau

    # loop over stages backwards
    for k = reverse(1:NS)
        # A
        z .= z .+ tab[:bᴵ, k].*Δt.*x
        y .= y .+ tab[:bᴱ, k].*Δt.*x
        # B
        v .= w # temporary
        sys(t + tab[:cᴱ, k]*Δt, stages[k], y, v)
        w .= w .+ v
        y .= 0
        # C
        z .= z .+ tab[:aᴵ, k, k].*Δt.*w
        y .= y .+ w
        w .= 0
        # D
        ImcAt!(sys, tab[:aᴵ, k, k]*Δt, z, v)
        s .= v .+ s
        z .= 0
        # E
        At_mul_B!(v, sys, s)
        y .= v .+ y
        s .= 0
        # F
        if k == 1
            x .= x .+ y
            y .= 0
        else
            x .= y .+ x
            z .= (tab[:aᴵ, k, k-1] .- tab[:bᴵ, k-1]).*Δt.*y .+ z
            y .= (tab[:aᴱ, k, k-1] .- tab[:bᴱ, k-1]).*Δt.*y
        end
    end

    return nothing
end