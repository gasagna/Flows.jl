# ---------------------------------------------------------------------------- #
# Type parameters:
# - X       : the state type
# - NS      : number of internal stages, used for dispatch with the stage cache
# - TYPETAG : :NORMAL or :LIN, use for dispatch too
# - ISADJ   : boolean for whether this is an adjoint integration
abstract type AbstractMethod{X, NS, TYPETAG, ISADJ} end

typetag(::AbstractMethod{X, NS, TYPETAG}) where {X, NS, TYPETAG} = TYPETAG

function _tag_map(tag::Symbol)
    tag == :NORMAL && return (:NORMAL,  false)
    tag == :TAN    && return (:LIN,     false)
    tag == :ADJ    && return (:LIN,      true)
    throw(ArgumentError("Type tag must be in (:NORMAL, :TAN, :ADJ)."))
end