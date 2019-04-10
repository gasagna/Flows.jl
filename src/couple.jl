using MacroTools: postwalk

export Coupled, couple, @all

# This is basically an N-tuple. We could use Base.Tuple{Any, Any}, obtaining 
# most of the functionality except `similar` and `copy`, which we would need to
# overload. That would be type piracy, so we create our implementation. The 
# only hacky bit is how to define the broadcast behaviour for the Coupled object.
# This used to be straightforward in julia v0.6, but we now need to use a macro 
# (@all, see below) to obtain the same behaviour. Since Coupled objects should only
# be used within this library code, this is not such a big issue.
struct Coupled{N, ARGS<:NTuple{N, Any}} <: AbstractVector{Float64}
    args::ARGS
    function Coupled(args::ARGS) where {N, ARGS<:NTuple{N, Any}}
        N < 100 || throw(ArgumentError("must be < 100"))
        return new{N, ARGS}(args)
    end
end

# Constructor
couple(args...) = Coupled(args)

# Extract parts using indexing (this is read-only)
Base.getindex(x::Coupled, i::Int) = x.args[i]
Base.length(x::Coupled) = length(x.args)

# Array interface
Base.similar(x::Coupled{N}) where {N} = couple(ntuple(i->similar(x[i]), N)...)
Base.copy(x::Coupled{N})    where {N} = couple(ntuple(i->copy(x[i]), N)...)

# Stuff for the @all macro
_isoperator(s::Symbol) = s ∈ (:+, :-, :*, :/, :.+, :.-, :.*, :./)
_isoperand(::Union{Expr, Number, QuoteNode}) = false
_isoperand(s::Symbol) = !_isoperator(s)

_getarg(x::Coupled, i::Int) = x[i]
_getarg(x,          i::Int) = x

# Broadcast any operations on all arguments of the Coupled object
macro all(expr)
    exprs = [postwalk(s->_isoperand(s) ? :(Flows._getarg($s, $i)) : s, expr) for i = 1:99]
    quote
        if typeof($(esc(expr.args[1]))) <: Coupled
            N = length($(esc(expr.args[1])))
            N >=  1 && $(esc(exprs[ 1]))
            N >=  2 && $(esc(exprs[ 2]))
            N >=  3 && $(esc(exprs[ 3]))
            N >=  4 && $(esc(exprs[ 4]))
            N >=  5 && $(esc(exprs[ 5]))
            N >=  6 && $(esc(exprs[ 6]))
            N >=  7 && $(esc(exprs[ 7]))
            N >=  8 && $(esc(exprs[ 8]))
            N >=  9 && $(esc(exprs[ 9]))
            N >= 10 && $(esc(exprs[10]))
            N >= 11 && $(esc(exprs[11]))
            N >= 12 && $(esc(exprs[12]))
            N >= 13 && $(esc(exprs[13]))
            N >= 14 && $(esc(exprs[14]))
            N >= 15 && $(esc(exprs[15]))
            N >= 16 && $(esc(exprs[16]))
            N >= 17 && $(esc(exprs[17]))
            N >= 18 && $(esc(exprs[18]))
            N >= 19 && $(esc(exprs[19]))
            N >= 20 && $(esc(exprs[20]))
            N >= 21 && $(esc(exprs[21]))
            N >= 22 && $(esc(exprs[22]))
            N >= 23 && $(esc(exprs[23]))
            N >= 24 && $(esc(exprs[24]))
            N >= 25 && $(esc(exprs[25]))
            N >= 26 && $(esc(exprs[26]))
            N >= 27 && $(esc(exprs[27]))
            N >= 28 && $(esc(exprs[28]))
            N >= 29 && $(esc(exprs[29]))
            N >= 30 && $(esc(exprs[30]))
            N >= 31 && $(esc(exprs[31]))
            N >= 32 && $(esc(exprs[32]))
            N >= 33 && $(esc(exprs[33]))
            N >= 34 && $(esc(exprs[34]))
            N >= 35 && $(esc(exprs[35]))
            N >= 36 && $(esc(exprs[36]))
            N >= 37 && $(esc(exprs[37]))
            N >= 38 && $(esc(exprs[38]))
            N >= 39 && $(esc(exprs[39]))
            N >= 40 && $(esc(exprs[40]))
            N >= 41 && $(esc(exprs[41]))
            N >= 42 && $(esc(exprs[42]))
            N >= 43 && $(esc(exprs[43]))
            N >= 44 && $(esc(exprs[44]))
            N >= 45 && $(esc(exprs[45]))
            N >= 46 && $(esc(exprs[46]))
            N >= 47 && $(esc(exprs[47]))
            N >= 48 && $(esc(exprs[48]))
            N >= 49 && $(esc(exprs[49]))
            N >= 50 && $(esc(exprs[50]))
            N >= 51 && $(esc(exprs[51]))
            N >= 52 && $(esc(exprs[52]))
            N >= 53 && $(esc(exprs[53]))
            N >= 54 && $(esc(exprs[54]))
            N >= 55 && $(esc(exprs[55]))
            N >= 56 && $(esc(exprs[56]))
            N >= 57 && $(esc(exprs[57]))
            N >= 58 && $(esc(exprs[58]))
            N >= 59 && $(esc(exprs[59]))
            N >= 60 && $(esc(exprs[60]))
            N >= 61 && $(esc(exprs[61]))
            N >= 62 && $(esc(exprs[62]))
            N >= 63 && $(esc(exprs[63]))
            N >= 64 && $(esc(exprs[64]))
            N >= 65 && $(esc(exprs[65]))
            N >= 66 && $(esc(exprs[66]))
            N >= 67 && $(esc(exprs[67]))
            N >= 68 && $(esc(exprs[68]))
            N >= 69 && $(esc(exprs[69]))
            N >= 70 && $(esc(exprs[70]))
            N >= 71 && $(esc(exprs[71]))
            N >= 72 && $(esc(exprs[72]))
            N >= 73 && $(esc(exprs[73]))
            N >= 74 && $(esc(exprs[74]))
            N >= 75 && $(esc(exprs[75]))
            N >= 76 && $(esc(exprs[76]))
            N >= 77 && $(esc(exprs[77]))
            N >= 78 && $(esc(exprs[78]))
            N >= 79 && $(esc(exprs[79]))
            N >= 80 && $(esc(exprs[80]))
            N >= 81 && $(esc(exprs[81]))
            N >= 82 && $(esc(exprs[82]))
            N >= 83 && $(esc(exprs[83]))
            N >= 84 && $(esc(exprs[84]))
            N >= 85 && $(esc(exprs[85]))
            N >= 86 && $(esc(exprs[86]))
            N >= 87 && $(esc(exprs[87]))
            N >= 88 && $(esc(exprs[88]))
            N >= 89 && $(esc(exprs[89]))
            N >= 90 && $(esc(exprs[90]))
            N >= 91 && $(esc(exprs[91]))
            N >= 92 && $(esc(exprs[92]))
            N >= 93 && $(esc(exprs[93]))
            N >= 94 && $(esc(exprs[94]))
            N >= 95 && $(esc(exprs[95]))
            N >= 96 && $(esc(exprs[96]))
            N >= 97 && $(esc(exprs[97]))
            N >= 98 && $(esc(exprs[98]))
            N == 99 && $(esc(exprs[99]))
        else typeof($(esc(expr.args[1])))
            $(esc(expr))
        end
        $(esc(expr.args[1]))
    end
end