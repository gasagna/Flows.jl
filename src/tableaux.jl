# ---------------------------------------------------------------------------- #
# Abstract tableau definition
abstract type AbstractTableau{T, NS} end

# number of stages
nstages(::AbstractTableau{T, NS}) where {T, NS} = NS

# ---------------------------------------------------------------------------- #
# Single tableau
struct Tableau{T, NS} <: AbstractTableau{T, NS}
    a::Matrix{T}
    b::Vector{T}
    e::Vector{T}
    c::Vector{T}
    function Tableau{NS}(a::Matrix{T},
                         b::Vector{T},
                         e::Vector{T},
                         c::Vector{T}) where {T, NS}
        NS == size(a, 1) == size(a, 2) == length(b) == length(c) == length(e) ||
            error("invalid tableau input")
        new{T, NS}(a, b, e, c)
    end
end

# outer constructor
function Tableau(a::Matrix, b::Vector, e::Vector, c::Vector)
    T = promote_type(eltype.((a, b, e, c))...)
    Tableau{length(c)}(convert(Matrix{T}, a), convert.(Vector{T}, (b, e, c))...)
end

# get coefficients of the tableau
Base.getindex(tab::Tableau, ::Symbol, i::Integer, j::Integer) = tab.a[i, j]
function Base.getindex(tab::Tableau{T}, t::Symbol, i::Integer)::T where {T}
    t == :b && return tab.b[i]
    t == :e && return tab.e[i]
    t == :c && return tab.c[i]
end

# Convert a tableau to have coefficient of given type
Base.convert(::Type{Tableau{T, NS}}, tab::Tableau{S, NS}) where {T, S, NS} =
    Tableau(convert(Matrix{T}, tab.a), convert(Vector{T}, tab.b),
            convert(Vector{T}, tab.e), convert(Vector{T}, tab.c))

# ~~ Tableau for IMEX schemes ~~~
struct IMEXTableau{T, NS} <: AbstractTableau{T, NS}
    tabI::Tableau{T, NS}
    tabE::Tableau{T, NS}
end

# outer constructor
function IMEXTableau(tabI::Tableau{TI, NS},
                     tabE::Tableau{TE, NS}) where {TI, TE, NS}
    T = promote_type(TI, TE)
    IMEXTableau(convert(Tableau{T, NS}, tabI), convert(Tableau{T, NS}, tabE))
end

Base.convert(::Type{IMEXTableau{T, NS}},
          tab::IMEXTableau{S, NS}) where {T, S, NS} =
    IMEXTableau(convert(Tableau{T, NS}, tab.tabI),
                convert(Tableau{T, NS}, tab.tabE))

# get coefficients of the tableau
function Base.getindex(tab::IMEXTableau{T},
                         t::Symbol,
                         i::Integer,
                         j::Integer)::T where {T}
    t == :aᴵ && return tab.tabI[:a, i, j]
    t == :aᴱ && return tab.tabE[:a, i, j]
    throw(ArgumentError("symbol $t not recognized"))
end

function Base.getindex(tab::IMEXTableau{T}, t::Symbol, i::Integer)::T where {T}
    t == :bᴵ && return tab.tabI[:b, i]
    t == :eᴵ && return tab.tabI[:e, i]
    t == :cᴵ && return tab.tabI[:c, i]
    t == :bᴱ && return tab.tabE[:b, i]
    t == :eᴱ && return tab.tabE[:e, i]
    t == :cᴱ && return tab.tabE[:c, i]
    throw(ArgumentError("symbol $t not recognized"))
end


# Tableaux from Cavaglieri and Bewley 2015

# ~ IMEXRKCB2
const CB2_I = Tableau([0//1  0//1  0//1;
                       0//1  2//5  0//1;
                       0//1  5//6  1//6],
                      [0//1, 5//6, 1//6],
                      [0//1, 4//5, 1//5],
                      [0//1, 2//5, 1//1])

const CB2_E = Tableau([0//1  0//1  0//1;
                       2//5  0//1  0//1;
                       0//1  1//1  0//1],
                      [0//1, 5//6, 1//6],
                      [0//1, 4//5, 1//5],
                      [0//1, 2//5, 1//1])

# ~ CB3e
const CB3e_I = Tableau([0//1  0//1   0//1  0//1;
                        0//1  1//3   0//1  0//1;
                        0//1  1//2   1//2  0//1;
                        0//1  3//4  -1//4  1//2],
                       [0//1, 3//4, -1//4, 1//2],
                       [0//1, 3//4, -1//4, 1//2],
                       [0//1, 1//3,  1//1, 1//1])

const CB3e_E = Tableau([0//1  0//1   0//1  0//1;
                        1//3  0//1   0//1  0//1;
                        0//1  1//1   0//1  0//1;
                        0//1  3//4   1//4  0//1],
                       [0//1, 3//4, -1//4, 1//2],
                       [0//1, 3//4, -1//4, 1//2],
                       [0//1, 1//3,  1//1, 1//1])

# ~ CB3c
const CB3c_I = Tableau([0//1   0//1                                              0//1                         0//1;
                        0//1   3375509829940//4525919076317                      0//1                         0//1;
                        0//1  -11712383888607531889907//32694570495602105556248  566138307881//912153721139   0//1;
                        0//1   673488652607//2334033219546                       493801219040//853653026979   184814777513//1389668723319],
                       [0//1,  673488652607//2334033219546,                      493801219040//853653026979,  184814777513//1389668723319],
                       [0//1,  366319659506//1093160237145,                      270096253287//480244073137,  104228367309//1017021570740],
                       [0//1,  3375509829940//4525919076317,                     272778623835//1039454778728, 1//1])

const CB3c_E = Tableau([0//1                          0//1                          0//1                           0//1;
                        3375509829940//4525919076317  0//1                          0//1                           0//1;
                        0//1                          272778623835//1039454778728   0//1                           0//1;
                        0//1                          673488652607//2334033219546   1660544566939//2334033219546   0//1],
                       [0//1,                         673488652607//2334033219546,  493801219040//853653026979,    184814777513//1389668723319],
                       [449556814708//1155810555193,  0//1,                         210901428686//1400818478499,   480175564215//1042748212601],
                       [0//1,                         3375509829940//4525919076317, 272778623835//1039454778728,   1//1])


# ~ CB4
const CB4_I = Tableau([0//1                          0//1                          0//1                          0//1                         0//1                        0//1;
                       1//8                          1//8                          0//1                          0//1                         0//1                        0//1;
                       216145252607//961230882893    257479850128//1143310606989   30481561667//101628412017     0//1                         0//1                        0//1;
                       232049084587//1377130630063  -381180097479//1276440792700  -54660926949//461115766612     344309628413//552073727558   0//1                        0//1;
                       232049084587//1377130630063   322009889509//2243393849156  -100836174740//861952129159   -250423827953//1283875864443  1//2                        0//1;
                       232049084587//1377130630063   322009889509//2243393849156  -195109672787//1233165545817  -340582416761//705418832319   463396075661//409972144477  323177943294//1626646580633],
                      [232049084587//1377130630063,  322009889509//2243393849156, -195109672787//1233165545817, -340582416761//705418832319,  463396075661//409972144477, 323177943294//1626646580633],
                      [5590918588//49191225249,      92380217342//122399335103,   -29257529014//55608238079,    -126677396901//66917692409,   384446411890//169364936833, 58325237543//207682037557],
                      [0,                            1//4,                         3//4,                         3//8,                        1//2,                       1//1])

const CB4_E = Tableau([0//1                          0//1                          0//1                          0//1                          0//1                         0//1;
                       1//4                          0//1                          0//1                          0//1                          0//1                         0//1;
                       153985248130//1004999853329   902825336800//1512825644809   0//1                          0//1                          0//1                         0//1;
                       232049084587//1377130630063   99316866929//820744730663     82888780751//969573940619     0//1                          0//1                         0//1;
                       232049084587//1377130630063   322009889509//2243393849156   57501241309//765040883867     76345938311//676824576433     0//1                         0//1;
                       232049084587//1377130630063   322009889509//2243393849156  -195109672787//1233165545817  -4099309936455//6310162971841  1395992540491//933264948679  0//1],
                      [232049084587//1377130630063,  322009889509//2243393849156, -195109672787//1233165545817, -340582416761//705418832319,   463396075661//409972144477,  323177943294//1626646580633],
                      [5590918588//49191225249,      92380217342//122399335103,   -29257529014//55608238079,    -126677396901//66917692409,    384446411890//169364936833,  58325237543//207682037557],
                      [0,                            1//4,                         3//4,                         3//8,                         1//2,                        1//1])

# PR to add more are welcome! Defaults to Float64: FIXME
const CB2  = convert(IMEXTableau{Float64, 3}, IMEXTableau(CB2_I,   CB2_E))
const CB3c = convert(IMEXTableau{Float64, 4}, IMEXTableau(CB3c_I, CB3c_E))
const CB3e = convert(IMEXTableau{Float64, 4}, IMEXTableau(CB3e_I, CB3e_E))
const CB4  = convert(IMEXTableau{Float64, 6}, IMEXTableau(CB4_I,   CB4_E))
