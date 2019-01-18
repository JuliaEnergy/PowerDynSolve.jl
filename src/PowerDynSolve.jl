# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

__precompile__()

module PowerDynSolve

using DifferentialEquations
using PowerDynBase
using Lazy: @>
using Parameters: @with_kw

# erros
include("Errors.jl")

# data types
include("GridProblems.jl")
include("GridSolutions.jl")

# methods
include("solve.jl")

export solve, operationpoint, GridProblem, tspan, TimeSeries
export realsolve, complexsolve

end # module PowerDynSolve
