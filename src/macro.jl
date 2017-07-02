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

function subst(expr::Expr, with)
    # Recursively scan an expression and if you find 
    # something like `foo[i]` replace with 
    # `foo.with[i]`, where with will be `x` or `q` 
    # for objects of AugmentedState type.
    if expr.head == :ref
        cur = expr.args[1]        # this is: :foo
        rep = Expr(:., cur, with)
        expr.args[1] = rep        # this is: :foo.with
    end
    for arg in expr.args
        subst(arg, with)
    end
    return expr
end
subst(expr, with) = nothing

# Broadcast expression to the `x` and `q` fields
# of the objects involved in the expression. This
# happens only if flag is true.
function broadcast2fields(flag::Bool, expr::Expr)
    if flag
        withx = subst(copy(expr), :(:x))
        withy = subst(copy(expr), :(:q))
        ret = quote 
                    $withx
                    $withy
              end
    else
        ret = expr
    end
    return ret
end