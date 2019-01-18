using Test

testlist = [
    ("intergration.jl", "Integration Tests"),
    ("gridsolutions.jl", "Grid Solutions Tests"),
    ("plotrecipes.jl", "Plot Recipes Tests"),
]

@testset "$desc" for (file, desc) in testlist
    @time include(file)
end
