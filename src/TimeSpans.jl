

abstract type AbstractTimeSpan end
abstract type AbstractSingleTimeSpan <: AbstractTimeSpan end
Base.convert(::Type{T}, ts) where {T <: AbstractSingleTimeSpan} = T(ts...)
Base.convert(::Type{T}, ts::T) where {T <: AbstractSingleTimeSpan} = ts

struct TimeSpan <: AbstractSingleTimeSpan
    tBegin::Time
    tEnd::Time
    function TimeSpan(tBegin::Time, tEnd::Time)
        @assert tBegin <= tEnd # We will always go forward in time.
        new(tBegin, tEnd)
    end
end
TimeSpan(tBegin, tEnd) = TimeSpan(convert(Time, tBegin), convert(Time, tEnd))

Base.in(t, tSpan::TimeSpan) = tSpan.tBegin <= t <= tSpan.tEnd
Base.convert(::Type{Tuple}, tSpan::TimeSpan) = (tSpan.tBegin, tSpan.tEnd)

abstract type AbstractMultipleTimeSpans <: AbstractTimeSpan end
# TODO: Implement Iteration Protocol for AbstractMultipleTimeSpans

struct MultipleTimeSpans <: AbstractMultipleTimeSpans
    tSpans::Tuple{Vararg{AbstractSingleTimeSpan}}
end

function areSortedTimeSpans(tss::Tuple{Vararg{AbstractSingleTimeSpan}}; warnOnly=false)
    if tss === () return true end # empty set of time spans is sorted
    t = tss[1].tBegin # the first time
    for (n, ts) in enumerate(tss)
        if t > ts.tBegin
            msg = "Incorrectly ordered or overlapping time spans: time span $(n-1) ends at $t and time span $n starts with $(ts.tBegin)"
            if warnOnly
                @warn msg
            else
                throw(GridSolutionError(msg))
            end
            return false
        end
        t = ts.tEnd
    end
    return true
end

struct OrderedTimeSpans <: AbstractMultipleTimeSpans
    tSpans::Tuple{Vararg{AbstractSingleTimeSpan}}
    function OrderedTimeSpans(tSpans::Tuple{Vararg{AbstractSingleTimeSpan}}; warnOnly=false)
        areSortedTimeSpans(tSpans, warnOnly=warnOnly)
        new(tSpans)
    end
end
OrderedTimeSpans(tss; kwargs...) = OrderedTimeSpans(map(ts -> convert(TimeSpan, ts), tss); kwargs...)

@inline Base.in(t, tSpans::AbstractMultipleTimeSpans) = any(map(ts -> in(t, ts), tSpans.tSpans))
@inline Base.getindex(tss::AbstractMultipleTimeSpans, i) = getindex(tss.tSpans, i)
@inline Base.lastindex(tss::AbstractMultipleTimeSpans) = lastindex(tss.tSpans)

@inline Base.convert(::Type{Tuple{Vararg{TimeSpan}}}, tSpans::AbstractMultipleTimeSpans) = tSpans.tSpans
@inline Base.convert(::Type{Tuple}, tSpans::AbstractMultipleTimeSpans) = map(tSpan -> convert(Tuple, tSpan), tSpans.tSpans)

@inline Base.convert(::Type{T}, tss) where {T <: AbstractMultipleTimeSpans} = T(map(ts -> convert(TimeSpan, ts), tss))

@inline Base.findfirst(t, tSpans::OrderedTimeSpans) = findfirst(ts -> t in ts, tSpans.tSpans)
