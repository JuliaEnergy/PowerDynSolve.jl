
abstract type AbstractMultipleGridSolutions <: AbstractGridSolution end

struct SubSolutionHandler
    sols::Tuple{Vararg{AbstractSingleGridSolutions}}
end
@inline solutionsOf(ssh::SubSolutionHandler) = ssh.sols

SubSolutionHandler(sols) = begin
    @show sols .|> typeof
    @show flattenGridSolutions(sols...) .|> typeof
    SubSolutionHandler(flattenGridSolutions(sols...))
end

# not performant but easy to write (:
flattenGridSolutions(sols::Vararg{AbstractSingleGridSolutions}) = sols
flattenGridSolutions(sol::AbstractSingleGridSolutions, sols::Vararg{AbstractGridSolution}) =
    (sol, flattenGridSolutions(sols...)...)
flattenGridSolutions(sol::AbstractMultipleGridSolutions, sols::Vararg{AbstractGridSolution}) =
    (solutionsOf(sol)..., flattenGridSolutions(sols...)...)

(s::SubSolutionHandler)(::Nothing, args...; missingIfNotFound = false, kwargs...) =
    if missingIfNotFound return missing
    else throw(GridSolutionError("Couldn't find entry in solution for args = $args."))
    end
(s::SubSolutionHandler)(n, args...; kwargs...) =
    s.sols[n](args...; kwargs...)

Base.convert(::Type{SubSolutionHandler}, sols::Tuple) = SubSolutionHandler(sols)

struct CompositeGridSolution <: AbstractMultipleGridSolutions
    ssh::SubSolutionHandler
    tSpans::OrderedTimeSpans
end
CompositeGridSolution(ssh::SubSolutionHandler) = CompositeGridSolution(ssh, convert(OrderedTimeSpans, map(tspan, ssh.sols)))
CompositeGridSolution(sols::Vararg{AbstractGridSolution}) = CompositeGridSolution(convert(SubSolutionHandler, sols))

# map over each element in the time array
(csol::CompositeGridSolution)(t::AbstractArray{Time}, args...; kwargs...) =
    csol.(t, args...; kwargs...)

# takes missingIfNotFound as a kw
# propagate the call to the correct subsolution
(csol::CompositeGridSolution)(t::Time, args...; kwargs...) = begin
    # find the solution (number) to which the time t belongs to
    solutionNumber = findfirst(t, timeSpansOf(csol))
    # call that solution via the SubSolutionHandler (to catch the case where no solutionNumber has been found)
    subSolutionHandlerOf(csol)(solutionNumber, t, args...; kwargs...)
end

@inline subSolutionHandlerOf(csol::CompositeGridSolution) = csol.ssh
@inline solutionsOf(csol::CompositeGridSolution) = csol |> subSolutionHandlerOf |> solutionsOf
@inline timeSpansOf(csol::CompositeGridSolution) = csol.tSpans
