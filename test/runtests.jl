using Test

testlist = [
    ("rootfunction.jl", "Root Function Tests"),
    ("operationpoint.jl", "Operation Point Tests"),
    ("intergration.jl", "Integration Tests"),
    ("gridsolutions.jl", "Grid Solutions Tests"),
    ("plotrecipes.jl", "Plot Recipes Tests"),
]

@testset "$desc" for (file, desc) in testlist
    @time include(file)
end
