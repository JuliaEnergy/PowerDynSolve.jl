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
U1 = 1+5im
S2 = -2+3im
y = 2+1.5im
parnodes = [SlackAlgebraic(U=U1), PQAlgebraic(S=S2)]
LY = [y -y; -y y ]
grid = GridDynamics(parnodes, LY)
start = State(grid, rand(SystemSize(grid)))
fp = @test_nowarn operationpoint(start)
root = PowerDynSolve.RootFunction(grid)
@test all(root(convert(AbstractVector{Float64}, fp)) .- zeros(SystemSize(grid)) .< 1e-8)
@test S2 ≈ fp[2, :s]
end

let
parnodes = [PQAlgebraic(S=2), PQAlgebraic(S=-2)]
LY = [im -im; -im im]
g = GridDynamics(parnodes, LY)
start = State(g, rand(SystemSize(g)))
fp = operationpoint(start)
end

let
parnodes = [SwingEq(H=2, P =2, D=1, Ω=50), PQAlgebraic(S=-2)]
LY = [im -im; -im im]
g = GridDynamics(parnodes, LY)
start = State(g, rand(SystemSize(g)))
fp = operationpoint(start)
end

let
U1 = 1+5im
P2 = -1
V2 = 2
parnodes = [SlackAlgebraic(U=U1), SwingEqLVS(H=1, P=P2, D=1, Ω=50, Γ=20, V=V2)]
LY = [im -im; -im im]
grid = GridDynamics(parnodes, LY)
start = State(grid, rand(SystemSize(grid)))
fp = @test_nowarn operationpoint(start)
root = PowerDynSolve.RootFunction(grid)
@test all(root(convert(AbstractVector{Float64}, fp)) .- zeros(SystemSize(grid)) .< 1e-8)
@test V2 ≈ fp[2, :v]
@test P2 ≈ fp[2, :p]
end
