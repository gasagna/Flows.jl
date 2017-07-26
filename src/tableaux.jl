import Base: getindex, heads, tails, convert

# abstract type
abstract type AbstractTableau{N, T} end

# number of stages
nstages(::AbstractTableau{N}) where {N} = N

# Tableau coefficients
struct Tableau{N, T} <: AbstractTableau{N, T}
    a::NTuple{N, NTuple{N, T}}  #
    b::NTuple{N, T}             #
    e::NTuple{N, T}             # embedded scheme
    c::NTuple{N, T}             #
    Tableau{N, T}(a::NTuple{N, NTuple{N, T}},
                  b::NTuple{N, T},
                  e::NTuple{N, T},
                  c::NTuple{N, T}) where {N, T} = new(a, b, e, c)
end

# Outer constructor: convert all coefficients to a common type
function Tableau(a::NTuple{N, NTuple{N, Number}},
                 b::NTuple{N, Number},
                 e::NTuple{N, Number},
                 c::NTuple{N, Number}) where N
    T = promote_tuple_type((a, b, e, c))
    Tableau{N, T}(convert_tuple(T, (a, b, e, c))...)
end

# Convert a tableau to have coefficient of given type
convert(::Type{Tableau{N, T}}, tab::Tableau{N, S}) where {N, T, S} =
    Tableau{N, T}(convert_tuple(T, (tab.a, tab.b, tab.e, tab.c))...)

# convert a tuple such that all elements have type T
convert_tuple(::Type{T}, a::Number) where T = T(a)
convert_tuple(::Type{T}, a::Tuple)  where T = (map(el->convert_tuple(T, el), a)...)

# find common type for promotion across elements of a tuple, maybe nested
 promote_tuple_type(h::Tuple)                 = _promote_tuple_type(h)
_promote_tuple_type(h::Tuple)                 = promote_type(_promote_tuple_type(heads(h)...), _promote_tuple_type(tails(h)...))
_promote_tuple_type(h::Tuple{Vararg{Number}}) = promote_type(map(typeof, h)...)

# get coefficients of the tableau
getindex(tab::Tableau,  ::Symbol, i::Integer, j::Integer) = tab.a[i][j]
function getindex(tab::Tableau, t::Symbol, i::Integer) 
    t == :b && return tab.b[i]
    t == :e && return tab.e[i]
    t == :c && return tab.c[i]
    throw(ArgumentError("symbol $t not recognized"))
end


# ~~ Tableau for IMEX schemes ~~~
struct IMEXTableau{N, T} <: AbstractTableau{N, T}
    tabᴵ::Tableau{N, T}
    tabᴱ::Tableau{N, T}
end

# outer constructor
function IMEXTableau(tabᴵ::Tableau{N, Tᴵ}, tabᴱ::Tableau{N, Tᴱ}) where {N, Tᴵ, Tᴱ} 
    T = promote_type(Tᴵ, Tᴱ)
    IMEXTableau{N, T}(convert(Tableau{N, T}, tabᴵ), convert(Tableau{N, T}, tabᴱ))
end

convert(::Type{IMEXTableau{N, T}}, tab::IMEXTableau{N, S}) where {N, T, S} =
    IMEXTableau{N, T}(map(t->convert(Tableau{N, T}, t), (tab.tabᴵ, tab.tabᴱ))...)

# get coefficients of the tableau
function getindex(tab::IMEXTableau,  t::Symbol, i::Integer, j::Integer)
    t == :aᴵ && return tab.tabᴵ[:a, i, j]
    t == :aᴱ && return tab.tabᴱ[:a, i, j]
    throw(ArgumentError("symbol $t not recognized"))
end
function getindex(tab::IMEXTableau, t::Symbol, i::Integer) 
    t == :bᴵ && return tab.tabᴵ[:b, i]
    t == :eᴵ && return tab.tabᴵ[:e, i]
    t == :cᴵ && return tab.tabᴵ[:c, i]
    t == :bᴱ && return tab.tabᴱ[:b, i]
    t == :eᴱ && return tab.tabᴱ[:e, i]
    t == :cᴱ && return tab.tabᴱ[:c, i]
    throw(ArgumentError("symbol $t not recognized"))
end

# Tableaux from Cavaglieri and Bewley 2015
export IMEXRKCB3c, IMEXRKCB3e, IMEXRKCB4

# ~ IMEXRKCB3e
const IMEXRKCB3e_I = Tableau(((0//1, 0//1,  0//1, 0//1),
                              (0//1, 1//3,  0//1, 0//1),
                              (0//1, 1//2,  1//2, 0//1),
                              (0//1, 3//4, -1//4, 1//2)),
                              (0//1, 3//4, -1//4, 1//2),
                              (0//1, 3//4, -1//4, 1//2),
                              (0//1, 1//3,  1//1, 1//1))

const IMEXRKCB3e_E = Tableau(((0//1, 0//1,  0//1, 0//1),
                              (1//3, 0//1,  0//1, 0//1),
                              (0//1, 1//1,  0//1, 0//1),
                              (0//1, 3//4,  1//4, 0//1)),
                              (0//1, 3//4, -1//4, 1//2),
                              (0//1, 3//4, -1//4, 1//2),
                              (0//1, 1//3,  1//1, 1//1))

const IMEXRKCB3e = IMEXTableau(IMEXRKCB3e_I, IMEXRKCB3e_E)

# ~ IMEXRKCB3c
const IMEXRKCB3c_I = Tableau(((0//1,  0//1,                                             0//1,                        0//1),
                              (0//1,  3375509829940//4525919076317,                     0//1,                        0//1),
                              (0//1, -11712383888607531889907//32694570495602105556248, 566138307881//912153721139,  0//1),
                              (0//1,  673488652607//2334033219546,                      493801219040//853653026979,  184814777513//1389668723319)),

                              (0//1,  673488652607//2334033219546,                      493801219040//853653026979,  184814777513//1389668723319),
                              (0//1,  366319659506//1093160237145,                      270096253287//480244073137,  104228367309//1017021570740),
                              (0//1,  3375509829940//4525919076317,                     272778623835//1039454778728, 1//1))

const IMEXRKCB3c_E = Tableau(((0//1,                         0//1,                         0//1,                          0//1),
                              (3375509829940//4525919076317, 0//1,                         0//1,                          0//1),
                              (0//1,                         272778623835//1039454778728,  0//1,                          0//1),
                              (0//1,                         673488652607//2334033219546,  1660544566939//2334033219546,  0//1)),

                              (0//1,                         673488652607//2334033219546,  493801219040//853653026979,    184814777513//1389668723319),
                              (449556814708//1155810555193,  0//1,                         210901428686//1400818478499,   480175564215//1042748212601),
                              (0//1,                         3375509829940//4525919076317, 272778623835//1039454778728,   1//1))

const IMEXRKCB3c = IMEXTableau(IMEXRKCB3c_I, IMEXRKCB3c_E)

# ~ IMEXRKCB4
const IMEXRKCB4_I = Tableau(((0//1,                         0//1,                         0//1,                         0//1,                        0//1,                       0//1), 
                             (1//8,                         1//8,                         0//1,                         0//1,                        0//1,                       0//1),   
                             (216145252607//961230882893,   257479850128//1143310606989,  30481561667//101628412017,    0//1,                        0//1,                       0//1),
                             (232049084587//1377130630063, -381180097479//1276440792700, -54660926949//461115766612,    344309628413//552073727558,  0//1,                       0//1),
                             (232049084587//1377130630063,  322009889509//2243393849156, -100836174740//861952129159,  -250423827953//1283875864443, 1//2,                       0//1),
                             (232049084587//1377130630063,  322009889509//2243393849156, -195109672787//1233165545817, -340582416761//705418832319,  463396075661//409972144477, 323177943294//1626646580633)),
                             (232049084587//1377130630063,  322009889509//2243393849156, -195109672787//1233165545817, -340582416761//705418832319,  463396075661//409972144477, 323177943294//1626646580633),
                             (5590918588//49191225249,      92380217342//122399335103,   -29257529014//55608238079,    -126677396901//66917692409,   384446411890//169364936833, 58325237543//207682037557),
                             (0,                            1//4,                         3//4,                         3//8,                        1//2,                       1//1))

const IMEXRKCB4_E = Tableau(((0//1,                         0//1,                         0//1,                         0//1,                         0//1,                        0//1), 
                             (1//4,                         0//1,                         0//1,                         0//1,                         0//1,                        0//1),   
                             (153985248130//1004999853329,  902825336800//1512825644809,  0//1,                         0//1,                         0//1,                        0//1),
                             (232049084587//1377130630063,  99316866929//820744730663,    82888780751//969573940619,    0//1,                         0//1,                        0//1),
                             (232049084587//1377130630063,  322009889509//2243393849156,  57501241309//765040883867,    76345938311//676824576433,    0//1,                        0//1),
                             (232049084587//1377130630063,  322009889509//2243393849156, -195109672787//1233165545817, -4099309936455//6310162971841, 1395992540491//933264948679, 0//1)),
                             (232049084587//1377130630063,  322009889509//2243393849156, -195109672787//1233165545817, -340582416761//705418832319,   463396075661//409972144477,  323177943294//1626646580633),
                             (5590918588//49191225249,      92380217342//122399335103,   -29257529014//55608238079,    -126677396901//66917692409,    384446411890//169364936833,  58325237543//207682037557),
                             (0,                            1//4,                         3//4,                         3//8,                         1//2,                        1//1))

const IMEXRKCB4 = IMEXTableau(IMEXRKCB4_I, IMEXRKCB4_E)

# PR to add more are welcome!