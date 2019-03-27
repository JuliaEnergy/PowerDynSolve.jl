begin
    using PowerDynBase
    using PowerDynSolve
    using Random
    using Test
    @show random_seed = 1234
    Random.seed!(random_seed)
end

################################################################################
# no indentation to easily run it in atom but still have the let environments
# for the actual test runs
################################################################################

# sub types of AbstractTimeSpan
let
_ts1 = (1, 2)
# ensure integers are used as input so the test for the type compatibility
# with DifferentialEquations.jl below makes sense
@test isa(_ts1, Tuple{Int, Int})

# TimeSpan
tspan1 = convert(PowerDynSolve.TimeSpan, _ts1)
# check types to be compatible with DifferentialEquations.jl
@test isa(tspan1.tBegin, PowerDynSolve.Time)
@test isa(tspan1.tEnd, PowerDynSolve.Time)
# check that only correct time spans are created
@test_throws MethodError convert(PowerDynSolve.TimeSpan, (1.0,))
@test_throws MethodError convert(PowerDynSolve.TimeSpan, 1.0)
@test_throws MethodError convert(PowerDynSolve.TimeSpan, (1.0, 2.0, 3.0))
# check methods
@test 1 ∈ tspan1 # left border
@test 2 ∈ tspan1 # right border
@test 1.5 ∈ tspan1 # within
@test 0.5 ∉ tspan1
@test 2.5 ∉ tspan1
@test convert(Tuple, tspan1) == (1., 2.)

# MultipleTimeSpans
tspan2 = convert(PowerDynSolve.MultipleTimeSpans, ((1, 2), (0.1, 1.5)))
# check methods
@test 1 ∈ tspan2
@test 2 ∈ tspan2
@test 1.5 ∈ tspan2
@test 0.5 ∈ tspan2 # is only in the second set
@test 0.1 ∈ tspan2
@test 0. ∉ tspan2
@test 2.5 ∉ tspan2
@test convert(Tuple, tspan2) == ((1., 2.), (0.1, 1.5))

# OrderedTimeSpans
# only ordered time spans are allowed
@test_throws PowerDynSolve.GridSolutionError convert(PowerDynSolve.OrderedTimeSpans, ((1, 2), (0.1, 1.5)))
#TODO: how to test:
# PowerDynSolve.OrderedTimeSpans(((1, 2), (0.1, 1.5)), warnOnly=true)
tspan3 = convert(PowerDynSolve.OrderedTimeSpans, ((1, 2), (2, 3), (4, 5)))
# check methods
@test 1 ∈ tspan3
@test 2 ∈ tspan3
@test 1.5 ∈ tspan3
@test 3 ∈ tspan3
@test 2.5 ∈ tspan3
@test 4 ∈ tspan3
@test 5 ∈ tspan3
@test 4.5 ∈ tspan3
@test 0.5 ∉ tspan3
@test 3.5 ∉ tspan3
@test 5.5 ∉ tspan3
@test convert(Tuple, tspan3) == ((1., 2.), (2., 3.), (4., 5))
@test tspan3[1] == PowerDynSolve.TimeSpan(1, 2)
@test tspan3[2] == PowerDynSolve.TimeSpan(2, 3)
@test tspan3[3] == PowerDynSolve.TimeSpan(4, 5)
@test tspan3[1:2] == (PowerDynSolve.TimeSpan(1, 2), PowerDynSolve.TimeSpan(2, 3))
@test tspan3[2:end] == (PowerDynSolve.TimeSpan(2, 3), PowerDynSolve.TimeSpan(4, 5))
@test tspan3[:] == (PowerDynSolve.TimeSpan(1, 2), PowerDynSolve.TimeSpan(2, 3), PowerDynSolve.TimeSpan(4, 5))
@test findfirst(1, tspan3) == 1
@test findfirst(1.5, tspan3) == 1
@test findfirst(2, tspan3) == 1
@test findfirst(2.5, tspan3) == 2
@test findfirst(3, tspan3) == 2
@test findfirst(4, tspan3) == 3
@test findfirst(4.5, tspan3) == 3
@test findfirst(5, tspan3) == 3
@test findfirst(0.5, tspan3) === nothing
@test findfirst(3.5, tspan3) === nothing
@test findfirst(5.5, tspan3) === nothing
end

# GridSolution
let
parnodes = [SwingEqLVS(H=1., P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1., P=-1, D=1, Ω=50, Γ=20, V=1)]
LY = [1.0im -im; -im im]
grid = GridDynamics(parnodes, LY)
start = State(grid, rand(SystemSize(grid)))
p = GridProblem(start, (0.,10.))
sol = solve(p)
@test (0., 10.) == tspan(sol)
@test (range(0, stop=10, length=1_000) .≈ tspan(sol, 1_000)) |> all
t = tspan(sol, 10_000)
for n=[1:2, :, 1, 2]
    for syms=[(:u,), (:p,), (:v,), (:int, 1), (:ω,) ]
        @test_nowarn sol(t, n, syms...)
    end
    @test sol(t, n, :int, 1) == sol(t, n, :ω)
end
@test_nowarn sol(t, :, :int, [1, 1]) # access the frequencies
end

# CompositeGridSolution of 2 GridSolution
let
function exampleRun(parnodes, ts)
    LY = [1.0im -im; -im im]
    grid = GridDynamics(parnodes, LY)
    start = State(grid, rand(SystemSize(grid)))
    p = GridProblem(start, ts)
    solve(p)
end

# TODO: test SubSolutionHandler
ts1 = (0.,10.)
ts2 = (15., 25.)
parnodes1 = [SwingEqLVS(H=1., P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1., P=-1, D=1, Ω=50, Γ=20, V=1)]
parnodes2 = [SwingEqLVS(H=1., P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1., P=0, D=1, Ω=50, Γ=20, V=1)]
sol1 = exampleRun(parnodes1, ts1)
sol2 = exampleRun(parnodes2, ts2)
ssh = PowerDynSolve.SubSolutionHandler((sol1, sol2))
@test_throws PowerDynSolve.GridSolutionError ssh(nothing)
@test ssh(nothing, missingIfNotFound=true) === missing

# TODO: test argument propagation to the solution of the ssh

csol = CompositeGridSolution(sol1, sol2)
@test_nowarn PowerDynSolve.SubSolutionHandler((csol, sol1))
@test_throws PowerDynSolve.GridSolutionError CompositeGridSolution(csol, sol1)
end
