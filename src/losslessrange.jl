# A lightweight range object that remembers the last element
struct LossLessRange{T, R<:AbstractRange{T}} <: AbstractRange{T}
        rng::R
       stop::T
    isLossy::Bool
    LossLessRange(rng::R, stop::T, isLossy::Bool) where {T, R} =
        new{T, R}(rng, stop, isLossy)
end

function LossLessRange(start, stop, step)
    # make sure we accept a positive step and that internally we 
    # change its sign depending on the order of start and stop
    step > 0 || throw(ArgumentError("Step must be positive. Got $step."))
    step = start > stop ? -step : step
    rng = start:step:stop
    last(rng) != stop ? LossLessRange(rng, oftype(first(rng), stop), true) :
                        LossLessRange(rng, oftype(first(rng), stop), false)
end

@inline Base.length(llr::LossLessRange) =
    llr.isLossy ? length(llr.rng) + 1 : length(llr.rng)

@inline Base.getindex(llr::LossLessRange, i::Integer) = 
    llr.isLossy ? (i == length(llr) ? llr.stop : llr.rng[i]) : llr.rng[i]