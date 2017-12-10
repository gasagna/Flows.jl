# A lightweight range object that remembers the last element
struct LossLessRange{T, R<:Range{T}} <: Range{T}
        rng::R
       stop::T
    isLossy::Bool
    LossLessRange(rng::R, stop::T, isLossy::Bool) where {T, R} =
        new{T, R}(rng, stop, isLossy)
end

function LossLessRange(start, stop, step)
    rng = start:step:stop
    last(rng) != stop ? LossLessRange(rng, oftype(first(rng), stop), true) :
                        LossLessRange(rng, oftype(first(rng), stop), false)
end

@inline Base.length(llr::LossLessRange) =
    llr.isLossy ? length(llr.rng) + 1 : length(llr.rng)

@inline Base.getindex(llr::LossLessRange, i::Integer) = 
    llr.isLossy ? (i == length(llr) ? llr.stop : llr.rng[i]) : llr.rng[i]