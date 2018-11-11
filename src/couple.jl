using MacroTools: postwalk

export Coupled, couple, @all

# This is basically an N-tuple. We could use Base.Tuple{Any, Any}, obtaining 
# most of the functionality except `similar` and `copy`, which we would need to
# overload. That would be type piracy, so we create our implementation.
struct Coupled{N, ARGS<:NTuple{N, Any}}
    args::ARGS
    function Coupled(args::ARGS) where {N, ARGS<:NTuple{N, Any}}
        N ∈ (1, 2, 3) || throw(ArgumentError("must be 1, 2 or 3 args"))
        return new{N, ARGS}(args)
    end
end

# constructor
couple(args...) = Coupled(args)

# extract parts (up to three)
Base.getindex(x::Coupled, i::Int) = x.args[i]
Base.length(x::Coupled) = length(x.args)

# array interface
Base.similar(x::Coupled{N}) where {N} = couple(ntuple(i->similar(x[i]), N)...)
Base.copy(x::Coupled{N})    where {N} = couple(ntuple(i->copy(x[i]), N)...)

# stuff for the macro
_isoperator(s::Symbol) = s ∈ (:+, :-, :*, :/, :.+, :.-, :.*, :./)
_isoperand(::Union{Expr, Number, QuoteNode}) = false
_isoperand(s::Symbol) = !_isoperator(s)

_getarg(x::Coupled, i::Int) = x[i]
_getarg(x,          i::Int) = x

# Broadcast any operations on all arguments of the Coupled object
macro all(expr)
    expr_1 = postwalk(s->_isoperand(s) ? :(Flows._getarg($s, 1)) : s, expr)
    expr_2 = postwalk(s->_isoperand(s) ? :(Flows._getarg($s, 2)) : s, expr)
    expr_3 = postwalk(s->_isoperand(s) ? :(Flows._getarg($s, 3)) : s, expr)
    quote
        if typeof($(esc(expr.args[1]))) <: Coupled{3}
            $(esc(expr_1))
            $(esc(expr_2))
            $(esc(expr_3))
        elseif typeof($(esc(expr.args[1]))) <: Coupled{2}
            $(esc(expr_1))
            $(esc(expr_2))
        else typeof($(esc(expr.args[1])))
            $(esc(expr))
        end
        $(esc(expr.args[1]))
    end
end