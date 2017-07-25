using MacroTools

export @over_i

# All of this is required because julia does not currently support using the dot 
# notation for broadcasting inside generated functions, which we use to generate
# the stepping code. See bug https://github.com/JuliaLang/julia/issues/21094.
# When this is solved (if it will be solved one day), all of this can be 
# eliminated and the dot notation can be used in the generated functions.
# The broadcasting behaviour for AugmentedState type will need to be defined too.#

macro over_i(expr)
    # Accepts expressions of the form
    # x[i] = ...
    @assert expr.head in [:(=), :(+=), :(-=), :(*=), :(/=)]
    @assert expr.args[1].head == :ref
    # the array referenced is the one on 
    # the left hand side, i.e. x in the example above
    arr = expr.args[1].args[1]
    # the index of the loop that need to vary
    ind = expr.args[1].args[2]
    quote
        @inbounds @simd for $(esc(ind)) in eachindex($(esc(arr)))
            $(esc(expr))
        end
    end
end

# Substitute expressions like `x[i]` with `x.with[i]`
function substitute(ex, with)
    MacroTools.postwalk(ex) do x
        return @capture(x, arr_[i]) ? :($arr.$with[i]) : x
    end
end

# Broadcast expression to the `x` and `q` fields of the
# `AugmentedState` objects involved in the expression `expr`.
function loop_i_sq(expr)
    quote 
        @over_i $(substitute(copy(expr), :x))
        @over_i $(substitute(copy(expr), :q))
    end
end

# Do not broadcast, and expect the types to have defined
# the linear indexing behaviour.
function loop_i_s(expr)
    quote
        @over_i $expr
    end
end