# ---------------------------------------------------------------------------- #
# Type parameters:
# - X       : the state type
# - NS      : number of internal stages, used for dispatch with the stage cache
# - TYPETAG : :NL or :LIN, use for dispatch too
# - ISADJ   : boolean
abstract type AbstractMethod{X, NS, TYPETAG, ISADJ} end


function _tag_map(tag::Symbol)
	tag == :NL  && return (:NL,  false)
	tag == :TAN && return (:LIN, false)
	tag == :ADJ && return (:LIN,  true)
	throw(ArgumentError("Type tag must be in (:NL, :TAN, :ADJ)."))
end