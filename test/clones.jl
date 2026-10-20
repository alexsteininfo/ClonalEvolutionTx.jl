@testset "clone table bookkeeping" begin
    ct = CloneTable()
    ClonalEvolutionTx._bind_rates!(ct, s -> 1.0 + s, s -> 0.5)
    c1 = ClonalEvolutionTx._add_clone!(ct, 0, Traits(0.0, 1.0), :founder, 0.0, 0)
    ClonalEvolutionTx._change_size!(ct, c1, 10)
    c2 = ClonalEvolutionTx._add_clone!(ct, c1, Traits(0.2, 0.3), :resistance, 1.0, 7)
    ClonalEvolutionTx._change_size!(ct, c2, 1)
    @test (c1, c2) == (1, 2) && n_clones(ct) == 2
    @test living_clones(ct) == [1, 2]
    @test ct.maxβ ≈ 1.2 && ct.maxδ == 0.5 && ct.maxα == 1.0
    @test ct.sumNβ ≈ 10 * 1.0 + 1.2 && ct.sumNδ ≈ 5.5 && ct.sumNα ≈ 10.3
    ClonalEvolutionTx._change_size!(ct, c1, -10)            # founder dies out
    @test living_clones(ct) == [2]
    @test ct.maxα == 0.3                                   # extinct clone left the bound
    @test ct.size[1] == 0 && n_clones(ct) == 2              # history kept
    @test traits(ct, 2) == Traits(0.2, 0.3) && ct.parent[2] == 1 && ct.kind[2] == :resistance
    @test_throws ErrorException ClonalEvolutionTx._change_size!(ct, c1, 1)   # no revival
end

@testset "binding rates" begin
    ct = CloneTable()
    ClonalEvolutionTx._bind_rates!(ct, s -> 1.0, s -> 0.5)
    ClonalEvolutionTx._add_clone!(ct, 0, Traits(-2.0, 1.0), :founder, 0.0, 0)
    @test_throws ArgumentError ClonalEvolutionTx._bind_rates!(ct, s -> 1.0 + s, s -> 0.5)   # β = -1
end

@testset "clone table edge cases" begin
    ct = CloneTable()                                      # default rates β = 1, δ = 0
    @test ct.maxβ == 0.0 && ct.maxδ == 0.0 && ct.maxα == 0.0 && isempty(living_clones(ct))
    a = ClonalEvolutionTx._add_clone!(ct, 0, Traits(0.0, 0.5), :founder, 0.0, 0)
    @test ct.β[a] == 1.0 && ct.δ[a] == 0.0 && ct.size[a] == 0
    @test ct.alive_pos[a] == 0
    @test_throws ErrorException ClonalEvolutionTx._change_size!(ct, a, -1)   # negative size
    @test ct.size[a] == 0
    # negative rate at clone creation
    ClonalEvolutionTx._bind_rates!(ct, s -> 1.0 + s, s -> 0.5)
    @test_throws ArgumentError ClonalEvolutionTx._add_clone!(ct, 0, Traits(-2.0, 1.0), :driver, 0.0, 1)
    @test n_clones(ct) == 1
    # swap-remove keeps positions consistent; maxima recomputed when holder dies
    ids = [ClonalEvolutionTx._add_clone!(ct, 0, Traits(0.1 * i, 0.1 * i), :driver, 0.0, i) for i in 1:4]
    for k in ids
        ClonalEvolutionTx._change_size!(ct, k, 3)
    end
    ClonalEvolutionTx._change_size!(ct, ids[1], -3)
    @test sort(living_clones(ct)) == ids[2:4]
    @test all(ct.alive_pos[k] == findfirst(==(k), ct.alive) for k in ct.alive)
    @test ct.alive_pos[ids[1]] == 0
    ClonalEvolutionTx._change_size!(ct, ids[4], -3)
    @test ct.maxβ ≈ 1.3 && ct.maxα ≈ 0.3
    @test ct.sumNβ ≈ 3 * (1.2 + 1.3) && ct.sumNα ≈ 3 * (0.2 + 0.3)
    @test living_clones(ct) isa Vector{Int32}
    # rebinding recomputes everything
    ClonalEvolutionTx._bind_rates!(ct, s -> 2.0, s -> 1.0)
    @test ct.maxβ == 2.0 && ct.sumNβ ≈ 12.0 && ct.sumNδ ≈ 6.0
end
