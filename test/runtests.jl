using Test

@time @testset "RootFunction" begin include("rootfunction.jl") end

@time @testset "Operation Point" begin include("operationpoint.jl") end

@time @testset "Integration" begin include("intergration.jl") end

@time @testset "Grid Solutions" begin include("gridsolutions.jl") end

@time @testset "Plot Recipes" begin include("plotrecipes.jl") end
