import MacroTools: @capture, postwalk, prettify

export codegen, step!

# simplify code, by removing operations that do nothing
function strip_step(expr::Expr)
    expr = postwalk(ex->@capture(ex, arg_ .+ 0.0*Δt_ .* z_) ? :($arg)   : ex, expr)
    expr = postwalk(ex->@capture(ex, arg_ .= arg_                ) ? quote end : ex, expr)
    return expr
end

# Generate code for a specific implementation
function codegen(impl::AbstractIMEXRKImplementation; strip::Bool=true)
    T       = typeof(impl)
    funcode = strip == true ? strip_step(_codegen(impl)) : _codegen(impl)
    quote
        function step!(scheme::IMEXRKScheme{S, $T}, 
                       g, 
                       A, 
                       t::Real, 
                       Δt::Real, 
                       x::S) where S
            $funcode
        end
    end
end

# generate methods for consistent combinations of implementations and tableaux
for embed in [true, false]
    for impl in [IMEXRK3R2R(IMEXRKCB3e, embed),
                 IMEXRK3R2R(IMEXRKCB3c, embed),
                 IMEXRK4R3R(IMEXRKCB4,  embed)]

        # create function
        code = codegen(impl) 
        @eval $code
    end
end