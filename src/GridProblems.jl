# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

using DiffEqBase: DEProblem, RECOMPILE_BY_DEFAULT
using PowerDynBase: AbstractState, OrdinaryGridDynamics, OrdinaryGridDynamicsWithMass
using LinearAlgebra

const iipfunc = true # is in-place function

struct GridProblem{P<:DEProblem, S<:AbstractState, T<:AbstractFloat} # T is for the timespan
    prob::P
    start::S
    timespan::Tuple{T, T}
end
function GridProblem(start::AbstractState{OrdinaryGridDynamics, V, T}, timespan; kwargs...) where {V,T}
    GridProblem(ODEProblem{iipfunc}(GridDynamics(start), convert(AbstractVector{V}, start), timespan; kwargs...),
        start, timespan)
end

function GridProblem(start::AbstractState{OrdinaryGridDynamicsWithMass, V, T}, timespan; kwargs...) where {V,T}
    odefunc = ODEFunction{iipfunc, RECOMPILE_BY_DEFAULT}(GridDynamics(start), mass_matrix=( start |> NetworkRHS |> masses .|> Int |> Diagonal ))
    GridProblem(ODEProblem{iipfunc}(odefunc, convert(AbstractVector{V}, start), timespan; kwargs...),
        start, timespan)
end

function GridProblem(g::G, start::AbstractState{G, V, T}, timespan; kwargs...) where {G, V, T}
    # G is the type of the grid dynamics, the function call already assures
    # type equality of the grid dynamics
    GridProblem(start, timespan; kwargs...)
end

function GridProblem(g::GridDynamics, start::AbstractVector, timespan; kwargs...)
    GridProblem(State(g, start), timespan)
end
