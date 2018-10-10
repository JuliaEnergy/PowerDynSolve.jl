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
