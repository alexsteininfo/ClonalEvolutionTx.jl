@testset "tree and cell mode: same cells for one seed" begin
    pt, pc = grown(:tree; N = 2000, seed = 7), grown(:cells; N = 2000, seed = 7)
    @test cell_ids(pt) == cell_ids(pc) && cell_clones(pt) == cell_clones(pc)
    @test mutations_per_cell(pt) == mutations_per_cell(pc)
    @test divisions_per_cell(pt) == divisions_per_cell(pc)
    @test pt.clones.size == pc.clones.size && pt.t == pc.t
end

@testset "neutral mutations do not change the dynamics" begin
    p0, p3 = grown(:tree; seed = 8, μ = 0.0), grown(:tree; seed = 8, μ = 3.0)
    @test cell_ids(p0) == cell_ids(p3) && p0.t == p3.t && p0.clones.size == p3.clones.size
    @test sum(mutations_per_cell(p0)) == 0 && sum(mutations_per_cell(p3)) > 0
end
