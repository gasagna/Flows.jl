export IMEXRKScheme

# Scheme definition. A scheme is defined by an implementation of a IMEXRK
# method and the typed storage required. The implementation is 
# part of the type parameters because it is used to dispatch to the 
# appropriate specialisation of the `step!` method. The internal 
# constructor allocates the required storage, based on the implementation
struct IMEXRKScheme{S, I<:AbstractIMEXRKImplementation}
    storage::Vector{S}
    function IMEXRKScheme{S, I}(x::S) where {S, I}
        new(S[similar(x) for i = 1:nstorage(I)])
    end
end

# Outer constructor. Might wish to augment the state for quadrature rules.
IMEXRKScheme(impl::AbstractIMEXRKImplementation, x, q...) =
    (z = aug_state(x, q...); IMEXRKScheme{typeof(z), typeof(impl)}(z))
