using MacroTools: postwalk

export Coupled, couple, @all, arg1, arg2

# This is basically a 2-tuple. We could use Base.Tuple{Any, Any}, obtaining 
# most of the functionality except `similar` and `copy`, which we would need to
# overload. That would be type piracy, so we create our implementation.
# 
# There are two use cases in this library for this type:
# 1) a type that wraps the main state vector `a::A` augmented
#    with an appropriate object `b::B` used for quadrature integration.
#    Note that the type is immutable, hence if a scalar function need to be
#    integrated, the quadrature value is stored in a one-element vector, 
#    in the field `b`
# 2) a type to integrate a couple system of equations, where the coupling
#    between the two parts is equivalent to a lower triangular coupling matrix.
#    This is useful for integrating the tangent equations jointly with the 
#    nonlinear equations.
struct Coupled{ARG1, ARG2}
    arg1::ARG1
    arg2::ARG2
end

# constructors
couple(arg1, arg2) = Coupled(arg1, arg2)

# extract parts
arg1(x::Coupled) = x.arg1
arg2(x::Coupled) = x.arg2
arg1(x::Any) = x
arg2(x::Any) = x

# stuff for the macro
_isoperator(s::Symbol) = s ∈ (:+, :-, :*, :/, :.+, :.-, :.*, :./)
_isoperand(::Union{Expr, Number, QuoteNode}) = false
_isoperand(s::Symbol) = !_isoperator(s)

macro all(expr)
    expr_1 = postwalk(s->_isoperand(s) ? :(arg1($s)) : s, expr)
    expr_2 = postwalk(s->_isoperand(s) ? :(arg2($s)) : s, expr)
    quote
        if typeof($(esc(expr.args[1]))) <: Coupled
            $(esc(expr_1))
            $(esc(expr_2))
        else
            $(esc(expr))
        end
        $(esc(expr.args[1]))
    end
end

Base.similar(x::Coupled) = couple(similar(arg1(x)), similar(arg2(x)))
Base.copy(   x::Coupled) = couple(   copy(arg1(x)),    copy(arg2(x)))