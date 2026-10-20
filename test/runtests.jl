using ClonalEvolutionTx, Test, Random, Statistics, Distributions, StableRNGs, AbstractTrees, TOML
include("fixtures.jl")

const ALL_TESTS = ["package", "traits", "clones", "doses", "tree", "population", "blocks",
                   "simulate", "invariants", "validation", "treatment", "schedules",
                   "measurements", "statistics", "sampling", "storage", "regression"]
const SELECTED = isempty(ARGS) ? filter(t -> isfile(joinpath(@__DIR__, t * ".jl")), ALL_TESTS) : ARGS

@testset "ClonalEvolutionTx.jl" begin
    for name in SELECTED
        @testset "$name" begin
            include(name * ".jl")
        end
    end
end
