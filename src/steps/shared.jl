# ---------------------------------------------------------------------------- #
# Type parameters:
# - X       : the state type
# - NS      : number of internal stages, used for dispatch with the stage cache
# - TYPETAG : :NORMAL or :LIN, use for dispatch too
# - ISADJ   : boolean for whether this is an adjoint integration
abstract type AbstractMethod{X, NS, TYPETAG, ISADJ} end

function typetag(::AbstractMethod{X, NS, TYPETAG, ISADJ}) where {X, NS, TYPETAG, ISADJ}
    TYPETAG == :NORMAL && return :NORMAL
    TYPETAG == :LIN    && ISADJ == false && return :TAN
    TYPETAG == :LIN    && ISADJ == true  && return :ADJ
end

function _tag_map(tag::Symbol)
    tag == :NORMAL && return (:NORMAL,  false)
    tag == :TAN    && return (:LIN,     false)
    tag == :ADJ    && return (:LIN,      true)
    throw(ArgumentError("Type tag must be in (:NORMAL, :TAN, :ADJ)."))
end