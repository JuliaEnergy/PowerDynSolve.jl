
abstract type AbstractMultipleGridSolutions <: AbstractGridSolution end

struct SubSolutionHandler
    sols::Tuple{Vararg{AbstractSingleGridSolutions}}
    SubSolutionHandler(sols) = new(flattenGridSolutions(sols...))
end
@inline solutionsOf(ssh::SubSolutionHandler) = ssh.sols

# not performant but easy to write (:
flattenGridSolutions(sols::Vararg{AbstractSingleGridSolutions}) = sols
flattenGridSolutions(sol::AbstractSingleGridSolutions, sols::Vararg{AbstractGridSolution}) =
    (sol, flattenGridSolutions(sols...)...)
flattenGridSolutions(sol::AbstractMultipleGridSolutions, sols::Vararg{AbstractGridSolution}) =
    (solutionsOf(sol)..., flattenGridSolutions(sols...)...)

function (s::SubSolutionHandler)(::Nothing, args...; missingIfNotFound::Bool=false, kwargs...)
    if missingIfNotFound
        missing
    else
        throw(GridSolutionError("Couldn't find entry in solution for args = $args."))
    end
end
# TODO: add issue to propagate missingIfNotFound to the Single Grid Solutions
(s::SubSolutionHandler)(num, args...; missingIfNotFound::Bool=false, kwargs...) =
    s.sols[num](args...; kwargs...)

Base.convert(::Type{SubSolutionHandler}, sols::Tuple) = SubSolutionHandler(sols)

struct CompositeGridSolution <: AbstractMultipleGridSolutions
    ssh::SubSolutionHandler
    tSpans::MultipleTimeSpans{true} # true => needs to be sorted
    function CompositeGridSolution(ssh::SubSolutionHandler, tSpans::MultipleTimeSpans)
        !issorted(tSpans) && throw(GridSolutionError("Please sort the grid solutions / time spans correctly."))
        new(ssh, tSpans)
    end
end
CompositeGridSolution(ssh::SubSolutionHandler) = begin
    CompositeGridSolution(ssh, convert(MultipleTimeSpans, map(tspan, ssh.sols)))
end
CompositeGridSolution(sols::Vararg{AbstractGridSolution}) =
    CompositeGridSolution(convert(SubSolutionHandler, sols))

# TODO: add issue on github
const ERROR_INDEXING_WITH_COLON = GridSolutionError("Indexing the nodes with `:` is not allowed for a CompositeGridSolution as there might be ambiguities. Please use something like 1:numberOfNodes instead.")
(csol::CompositeGridSolution)(t, ::Colon, args...; kwargs...) =
    throw(ERROR_INDEXING_WITH_COLON)

# map over each element in the time array
(csol::CompositeGridSolution)(t, n, args...; kwargs...) = begin
    csol.(t', n, args...; kwargs...) # note that the ' makes it an outer/product iteration
end
(csol::CompositeGridSolution)(t, n::Number, args...; kwargs...) = begin
    csol.(t, n, args...; kwargs...) # note that the ' makes it an outer/product iteration
end

# takes missingIfNotFound as a kw
# propagate the call to the correct subsolution
(csol::CompositeGridSolution)(t::Number, n::Number, args...; kwargs...) = begin
    # @show t, args
    # find the solution (number) to which the time t belongs to
    solutionNumber = findfirst(t, timeSpansOf(csol))
    # @show solutionNumber
    # call that solution via the SubSolutionHandler (to catch the case where no solutionNumber has been found)
    subSolutionHandlerOf(csol)(solutionNumber, t, n, args...; kwargs...)
end

@inline subSolutionHandlerOf(csol::CompositeGridSolution) = csol.ssh
@inline solutionsOf(csol::CompositeGridSolution) = csol |> subSolutionHandlerOf |> solutionsOf
@inline timeSpansOf(csol::CompositeGridSolution) = csol.tSpans
