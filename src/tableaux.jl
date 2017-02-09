export Tableau, 
       IMEXTableau, 
       nstages

# TODO
# ~ add isembedded

# Tableau coefficients are defined in the type parameter
immutable Tableau{a_, b_, b̂_, c_} 
    # N enforces consistent sizes of the input
    Tableau{N}(a::NTuple{N, NTuple{N}}, b::NTuple{N}, b̂::NTuple{N}, c::NTuple{N}) = new()
end

Tableau(a, b, b̂, c) = Tableau{a, b, b̂, c}(a, b, b̂, c)

# number of stages
nstages{a, b, b̂, c}(::Type{Tableau{a, b, b̂, c}}) = length(a)

# get coefficients of the tableau
Base.getindex{a, b, b̂, c}(::Type{Tableau{a, b, b̂, c}}, ::Type{Val{:a}}, i::Integer, j::Integer) = a[i][j]
Base.getindex{a, b, b̂, c}(::Type{Tableau{a, b, b̂, c}}, ::Type{Val{:b}}, i::Integer)             = b[i]
Base.getindex{a, b, b̂, c}(::Type{Tableau{a, b, b̂, c}}, ::Type{Val{:b̂}}, i::Integer)             = b̂[i]
Base.getindex{a, b, b̂, c}(::Type{Tableau{a, b, b̂, c}}, ::Type{Val{:c}}, i::Integer)             = c[i]


# Tableau for IMEX scheme
immutable IMEXTableau{Tᴵ<:Tableau, Tᴱ<:Tableau} end

IMEXTableau{T<:Tableau, S<:Tableau}(tᴵ::T, tᴱ::S) = IMEXTableau{typeof(tᴵ), typeof(tᴱ)}()

# get coefficients of the tableau
Base.getindex{Tᴵ<:Tableau, Tᴱ<:Tableau}(t::Type{IMEXTableau{Tᴵ, Tᴱ}}, ::Type{Val{:aᴵ}}, i::Integer, j::Integer) = Tᴵ[Val{:a}, i, j]
Base.getindex{Tᴵ<:Tableau, Tᴱ<:Tableau}(t::Type{IMEXTableau{Tᴵ, Tᴱ}}, ::Type{Val{:bᴵ}}, i::Integer) = Tᴵ[Val{:b}, i]
Base.getindex{Tᴵ<:Tableau, Tᴱ<:Tableau}(t::Type{IMEXTableau{Tᴵ, Tᴱ}}, ::Type{Val{:b̂ᴵ}}, i::Integer) = Tᴵ[Val{:b̂}, i]
Base.getindex{Tᴵ<:Tableau, Tᴱ<:Tableau}(t::Type{IMEXTableau{Tᴵ, Tᴱ}}, ::Type{Val{:cᴵ}}, i::Integer) = Tᴵ[Val{:c}, i]
Base.getindex{Tᴵ<:Tableau, Tᴱ<:Tableau}(t::Type{IMEXTableau{Tᴵ, Tᴱ}}, ::Type{Val{:aᴱ}}, i::Integer, j::Integer) = Tᴱ[Val{:a}, i, j]
Base.getindex{Tᴵ<:Tableau, Tᴱ<:Tableau}(t::Type{IMEXTableau{Tᴵ, Tᴱ}}, ::Type{Val{:bᴱ}}, i::Integer) = Tᴱ[Val{:b}, i]
Base.getindex{Tᴵ<:Tableau, Tᴱ<:Tableau}(t::Type{IMEXTableau{Tᴵ, Tᴱ}}, ::Type{Val{:b̂ᴱ}}, i::Integer) = Tᴱ[Val{:b̂}, i]
Base.getindex{Tᴵ<:Tableau, Tᴱ<:Tableau}(t::Type{IMEXTableau{Tᴵ, Tᴱ}}, ::Type{Val{:cᴱ}}, i::Integer) = Tᴱ[Val{:c}, i]