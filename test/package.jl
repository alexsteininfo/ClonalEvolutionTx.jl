@testset "package identity" begin
    proj = TOML.parsefile(joinpath(pkgdir(ClonalEvolutionTx), "Project.toml"))
    @test proj["name"] == "ClonalEvolutionTx"
    @test proj["uuid"] != "7b855ee6-6887-412f-a571-26d20a5a92d7"
    @test Set(keys(proj["deps"])) ⊆ Set(["AbstractTrees", "Distributions", "Random", "StableRNGs", "Statistics"])
end
