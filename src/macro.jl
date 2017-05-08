export @over_i

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