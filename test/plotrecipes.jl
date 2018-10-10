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

let
# test the helper functions
@test PowerDynSolve.tslabel(:bla, 3) == "bla3"
@test PowerDynSolve.tslabel(:bla, [3,4]) == ["bla3", "bla4"]
@test PowerDynSolve.tslabel(:bla, 3, 4) == "bla3_4"
@test PowerDynSolve.tslabel(:bla, [3, 2], 4) == ["bla3_4", "bla2_4"]
@test PowerDynSolve.tstransform(zeros(4)) == zeros(4)
@test PowerDynSolve.tstransform(zeros(2, 6)) == zeros(6, 2)
@test_throws MethodError PowerDynSolve.tstransform(zeros(2,4,6))
end

let
using RecipesBase
KW = Dict{Symbol, Any}
RecipesBase.is_key_supported(::Symbol) = true
# generate some solution
parnodes = [SwingEqLVS(H=1., P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1., P=-1, D=1, Ω=50, Γ=20, V=1)]
LY = [1.0im -im; -im im]
grid = GridDynamics(parnodes, LY)
start = State(grid, rand(SystemSize(grid)))
p = GridProblem(start, (0.,10.))
sol = solve(p)
for callargs in [(:v,), (:φ,), (:p,), (:q,), (:int, 1), (:ω,)]
    data_list = RecipesBase.apply_recipe(KW(), sol, : , callargs...)
    @test length(data_list) == 1
    plotargs = data_list[1].args
    @test length(plotargs) == 2 + length(callargs)
    @test plotargs[1] === sol
    @test plotargs[2] == Base.OneTo(length(parnodes)) # range of all nodes
    @test plotargs[3:end] == callargs
    data_list = RecipesBase.apply_recipe(KW(), plotargs...)
    @test length(data_list) == 1
    plotargs = data_list[1].args
    @test plotargs[1] == tspan(sol, PowerDynSolve.PLOT_TTIME_RESOLUTION)
    @test plotargs[2] == PowerDynSolve.tstransform(sol(tspan(sol, PowerDynSolve.PLOT_TTIME_RESOLUTION), 1:length(parnodes), callargs...))
    @test size(plotargs[2]) == (PowerDynSolve.PLOT_TTIME_RESOLUTION, length(parnodes))

    nodenum = 1
    data_list = RecipesBase.apply_recipe(KW(), sol, nodenum, callargs...)
    @test length(data_list) == 1
    plotargs = data_list[1].args
    @test plotargs[1] == tspan(sol, PowerDynSolve.PLOT_TTIME_RESOLUTION)
    plotargs[2]
    sol(tspan(sol, PowerDynSolve.PLOT_TTIME_RESOLUTION), nodenum, callargs...)
    @test all(plotargs[2] .≈ sol(tspan(sol, PowerDynSolve.PLOT_TTIME_RESOLUTION), nodenum, callargs...))
end
end
