# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

using PowerDynBase: AbstractState, GridDynamics, State
using NLsolve

"""
    struct RootFunction
        grid::GridDynamics
    end

Basic data structure providing a method to evalute a subtype of [`PowerDynSolve.GridDynamics`](@ref) as
a right-hand-side function that can be used for root searching.
"""
struct RootFunction
    grid::GridDynamics
end
"""
    function (r::RootFunction)(x_out, x_in)

Evaluate the power grid dynamics of `r` for `x_in` and write it in `x_out`.
"""
function (r::RootFunction)(x_out, x_in)
    r.grid(x_out, x_in, 0., 0.)
end
"""
    function (r::RootFunction)(x_in)

Evaluate the power grid dynamics of `r` for `x_in` and return the result.
"""
function (r::RootFunction)(x_in::AbstractVector)
    x_out = similar(x_in)
    r.grid(x_out, x_in, 0., 0.)
    x_out
end
function (r::RootFunction)(s::AbstractState{G, V, T}) where {G, V, T}
    r(convert(Vector{V}, s))
end

"""
    function operationpoint(start::AbstractState)
    function operationpoint(grid::GridDynamics, start::AbstractState)
    function operationpoint(grid::GridDynamics, start::AbstractVector)

Find the operation point (fixed point of the right-hand-side) of a power grid (either `grid` or inferred)
by using the initial value `start`.
"""
function operationpoint(start::AbstractState{G, V, T}) where {G, V, T}
    if SlackAlgebraic ∉ start |> Nodes .|> parametersof .|> typeof
        @warn "currently not making any checks concerning assumptions
                of whether its possible to find the fixed point"
    end
    if SwingEq ∈ start |> Nodes .|> parametersof .|> typeof
        @warn "found SwingEq node but these should be SwingEqLVS (just SwingEq is not yet supported for operation point search)"
    end
    grid = GridDynamics(start)
    rootfct = RootFunction(grid)
    initp = convert(AbstractVector{V}, start)
    res = nlsolve(rootfct, initp)
    if ~res.f_converged
        throw(OperationPointError("solver did not converge"))
    end
    State(grid, res.zero)
end
function operationpoint(grid::GridDynamics, start::AbstractState)
    @assert grid === GridDynamics(start)
    operationpoint(start)
end
function operationpoint(grid::GridDynamics, start::AbstractVector)
    operationpoint(State(grid, start))
end
