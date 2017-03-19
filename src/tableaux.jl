export Tableau, 
       IMEXTableau, 
       nstages,
       IMEXRKCB3c, 
       IMEXRKCB3e 

import Base: convert, getindex, eltype

# TODO
# ~ add isembedded

# Tableau coefficients are defined in the type parameter
struct Tableau{T, a, b, b̂, c, N} 
    Tableau{T, a, b, b̂, c, N}(::NTuple{N, NTuple{N, T}}, 
                              ::NTuple{N, T}, 
                              ::NTuple{N, T}, 
                              ::NTuple{N, T}) where {T, a, b, b̂, c, N} = new()
end

# find common type that holds all elements
promote_tuple_type(a::Number) = typeof(a)
promote_tuple_type(a::Tuple) = promote_type(map(promote_tuple_type, a)...)

# convert a tuple such that all elements have type T
convert_tuple{T}(::Type{T}, a::Number) = T(a)
convert_tuple{T}(::Type{T}, a::Tuple) = (map(el->convert_tuple(T, el), a)...)

function Tableau{N}(a::NTuple{N, NTuple{N, Number}}, 
                    b::NTuple{N, Number}, 
                    b̂::NTuple{N, Number}, 
                    c::NTuple{N, Number})
    T = promote_tuple_type((a, b, b̂, c))
    Tableau{T, a, b, b̂, c, N}(convert_tuple(T, (a, b, b̂, c))...)
end

function convert{S, T, a, b, b̂, c, N}(::Type{Tableau{S}}, ::Tableau{T, a, b, b̂, c, N}) 
    tab = convert_tuple(S, (a, b, b̂, c))
    Tableau{S, tab..., N}(tab...)
end

# 
eltype{T, a, b, b̂, c, N}(::Type{Tableau{T, a, b, b̂, c, N}}) = T

# number of stages
nstages{T, a, b, b̂, c, N}(::Type{Tableau{T, a, b, b̂, c, N}}) = N

# get coefficients of the tableau
getindex{T, a, b, b̂, c, N}(::Type{Tableau{T, a, b, b̂, c, N}}, ::Type{Val{:a}}, i::Integer, j::Integer) = a[i][j]
getindex{T, a, b, b̂, c, N}(::Type{Tableau{T, a, b, b̂, c, N}}, ::Type{Val{:b}}, i::Integer)             = b[i]
getindex{T, a, b, b̂, c, N}(::Type{Tableau{T, a, b, b̂, c, N}}, ::Type{Val{:b̂}}, i::Integer)             = b̂[i]
getindex{T, a, b, b̂, c, N}(::Type{Tableau{T, a, b, b̂, c, N}}, ::Type{Val{:c}}, i::Integer)             = c[i]

# Tableau for IMEX scheme
immutable IMEXTableau{Tᴵ<:Tableau, Tᴱ<:Tableau} end
IMEXTableau{T<:Tableau, S<:Tableau}(tᴵ::T, tᴱ::S) = IMEXTableau{typeof(tᴵ), typeof(tᴱ)}()

# number of stages
nstages{Tᴵ<:Tableau, Tᴱ<:Tableau}(::Type{IMEXTableau{Tᴵ, Tᴱ}}) = nstages(Tᴵ)

# get coefficients of the tableau
getindex{Tᴵ<:Tableau, Tᴱ<:Tableau}(t::Type{IMEXTableau{Tᴵ, Tᴱ}}, ::Type{Val{:aᴵ}}, i::Integer, j::Integer) = Tᴵ[Val{:a}, i, j]
getindex{Tᴵ<:Tableau, Tᴱ<:Tableau}(t::Type{IMEXTableau{Tᴵ, Tᴱ}}, ::Type{Val{:bᴵ}}, i::Integer) = Tᴵ[Val{:b}, i]
getindex{Tᴵ<:Tableau, Tᴱ<:Tableau}(t::Type{IMEXTableau{Tᴵ, Tᴱ}}, ::Type{Val{:b̂ᴵ}}, i::Integer) = Tᴵ[Val{:b̂}, i]
getindex{Tᴵ<:Tableau, Tᴱ<:Tableau}(t::Type{IMEXTableau{Tᴵ, Tᴱ}}, ::Type{Val{:cᴵ}}, i::Integer) = Tᴵ[Val{:c}, i]
getindex{Tᴵ<:Tableau, Tᴱ<:Tableau}(t::Type{IMEXTableau{Tᴵ, Tᴱ}}, ::Type{Val{:aᴱ}}, i::Integer, j::Integer) = Tᴱ[Val{:a}, i, j]
getindex{Tᴵ<:Tableau, Tᴱ<:Tableau}(t::Type{IMEXTableau{Tᴵ, Tᴱ}}, ::Type{Val{:bᴱ}}, i::Integer) = Tᴱ[Val{:b}, i]
getindex{Tᴵ<:Tableau, Tᴱ<:Tableau}(t::Type{IMEXTableau{Tᴵ, Tᴱ}}, ::Type{Val{:b̂ᴱ}}, i::Integer) = Tᴱ[Val{:b̂}, i]
getindex{Tᴵ<:Tableau, Tᴱ<:Tableau}(t::Type{IMEXTableau{Tᴵ, Tᴱ}}, ::Type{Val{:cᴱ}}, i::Integer) = Tᴱ[Val{:c}, i]


# Tableaux from Cavaglieri and Bewley 2015
# ~ RKW3
const RKW3 = Tableau(((0//1,  0//1,  0//1),
                      (8//15, 0//1 , 0//1),
                      (1//4,  5//12, 0//1)), 
                      (1//4,  0//1,  3//4),
                      (1//4,  0//1,  3//4),
                      (0//1,  8//15, 2//3))

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