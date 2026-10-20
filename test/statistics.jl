@testset "SFS and burdens on a hand-built tree" begin
    pop = small_tree()
    @test site_frequency_spectrum(pop) == [4, 3]       # k=1: 1+3; k=2: 2 + 1 clonal
    @test mutations_per_cell(pop) == [4, 6] && total_mutational_burden(pop) == 7
    @test clonal_mutations(pop) == 3 && mean_burden(pop) == 5.0
    @test sum(k * s for (k, s) in enumerate(site_frequency_spectrum(pop))) == sum(mutations_per_cell(pop))
    @test_throws ArgumentError site_frequency_spectrum(initialize_population(2; mode = :cells))
end

@testset "expected SFS of exponential growth (Stein & Werner 2025, eq. 9)" begin
    for (d, runs) in ((0.0, 40), (0.5, 40))
        block = BirthDeathBlock(birth = s -> 1.0, death = s -> d, μ = 1.0, stop_above = 10_000, restart_on_extinction = true)
        sfs = mean(site_frequency_spectrum(simulate!(initialize_population(1), block, StableRNG(s)))[1:5] for s in 1:runs)
        expected = [2 * 10_000 * sum((d / 1.0)^l / ((k + l) * (k + l + 1)) for l in 0:200) for k in 1:5]
        @test all(isapprox.(sfs, expected; rtol = 0.05))
    end
end

@testset "clone sizes and nested mutation frequencies" begin
    pop = initialize_population(5; mode = :cells)
    ct = pop.clones          # clone-level statistics read only the table, so no cells are moved
    A, C = ClonalEvolutionTx._add_clone!, ClonalEvolutionTx._change_size!
    c2 = A(ct, 1, Traits(0.0, 0.5), :resistance, 1.0, 1);  C(ct, c2, 3)
    c3 = A(ct, c2, Traits(0.0, 0.4), :resistance, 2.0, 2); C(ct, c3, 2)
    c4 = A(ct, 1, Traits(0.1, 1.0), :driver, 3.0, 3);      C(ct, c4, 4)
    @test clone_sizes(pop) == [5, 3, 2, 4]
    @test mutation_frequencies(pop) == [5, 2, 4]        # clone 2 includes its child clone 3
end

@testset "repainting neutral mutations" begin
    pop = grown(:tree; N = 2000, μ = 1.0)
    ids = copy(cell_ids(pop))
    repaint_neutral!(pop, 0.0, StableRNG(1))
    @test all(==(0), mutations_per_cell(pop))
    # one repaint has ~6% relative noise (early draws are shared by many cells): average 25
    ratios = map(2:26) do s
        repaint_neutral!(pop, 2.0, StableRNG(s))
        mean(mutations_per_cell(pop)) / (2 * mean(divisions_per_cell(pop)))
    end
    @test isapprox(mean(ratios), 1; rtol = 0.05)
    @test cell_ids(pop) == ids
    @test_throws ArgumentError repaint_neutral!(grown(:cells), 1.0, StableRNG(3))
end

@testset "Hill numbers" begin
    c = [5, 3, 2]
    @test hill_number(c, 0) == 3
    @test hill_number(c, 2) ≈ 1 / 0.38
    @test hill_number(c, Inf) ≈ 2.0
    @test hill_number(c, 1) ≈ exp(-sum(p * log(p) for p in (0.5, 0.3, 0.2)))
    @test hill_number(c, 1 + 1e-8) ≈ hill_number(c, 1) rtol = 1e-6
    @test renyi_entropy(c, 2) ≈ log(1 / 0.38)
end

@testset "Hill numbers from an SFS: both definitions" begin
    sfs = [4, 0, 3]                     # N = 3, A = 4·1 + 3·3 = 13, mean burden 13/3
    @test sfs_hill_number(sfs, 0) == 7 == sfs_hill_number(sfs, 0; frequencies = :cell_fraction)
    @test sfs_hill_number(sfs, 2) ≈ 169 / 31
    @test sfs_hill_number(sfs, Inf) ≈ 13 / 3
    @test sfs_hill_number(sfs, 2; frequencies = :cell_fraction) ≈ 9 / 31
    @test sfs_renyi_entropy(sfs, 1; frequencies = :cell_fraction) ≈ (4 / 3) * log(3)
    @test sfs_hill_number(sfs, Inf; frequencies = :cell_fraction) ≈ 1.0
    m̄ = 13 / 3      # conversions of SPEC §8
    @test sfs_hill_number(sfs, 2) ≈ m̄^2 * sfs_hill_number(sfs, 2; frequencies = :cell_fraction)
    @test sfs_renyi_entropy(sfs, 1) ≈ sfs_renyi_entropy(sfs, 1; frequencies = :cell_fraction) / m̄ + log(m̄)
    @test_throws ArgumentError sfs_hill_number(sfs, 2; frequencies = :foo)
end

@testset "statistic types inside measurements" begin
    pop = grown(:tree; N = 300)
    r = measure(pop, [SFS(), CloneSizes(), HillNumbers(), SFSHillNumbers(q = [0, 2]), SensitiveResistant()])
    @test r.sfs == site_frequency_spectrum(pop) && r.clone_sizes == clone_sizes(pop)
    @test r.sfs_hill_numbers[1] == total_mutational_burden(pop)
    @test sum(r.sensitive_resistant) == popsize(pop)
end

@testset "edge cases and distances" begin
    @test hill_number([0, 0], 2) == 0.0 && renyi_entropy([0, 0], 1) == -Inf
    @test sfs_hill_number([0, 0, 0], 1) == 0.0
    @test_throws ArgumentError hill_number([1, 2], NaN)
    @test hill_number([1, 3], -1) ≈ (0.25^-1 + 0.75^-1)^(1 / 2)
    pop = small_tree()
    @test pairwise_distances(pop) == [4]                # 1 + 3 below the MRCA
    @test coalescence_times(pop) == [pop.t - 1.0]
    @test sfs_hill_number(Int[], 2) == 0.0
    empty = initialize_population(1)
    @test_throws ArgumentError pairwise_distances(initialize_population(2; mode = :cells))
    @test_throws ArgumentError repaint_neutral!(small_tree(), -1.0, StableRNG(1))
    @test n_sensitive(pop) + n_resistant(pop) == popsize(pop)
    cellpop = grown(:cells; N = 200)
    @test sum(clone_sizes(cellpop)) == popsize(cellpop) && mean_burden(cellpop) > 0
    @test measure(MutationsPerCell(), pop) == [4, 6] && statistic_name(Divisions()) == :divisions
    big = grown(:tree; N = 60)
    @test length(pairwise_distances(big)) == 60 * 59 ÷ 2
    @test all(>=(0), coalescence_times(big))
    @test clonal_mutations(big) <= minimum(mutations_per_cell(big))
end

@testset "Hill numbers at q = -Inf and large |q|" begin   # final review M1
    c = [5, 3, 2]                                     # p = 0.5, 0.3, 0.2
    @test hill_number(c, -Inf) ≈ 5.0                  # 1 / min p
    @test isapprox(hill_number(c, -1e6), 5.0; rtol = 1e-4) && isapprox(hill_number(c, -100), 5.0; rtol = 0.02)
    @test isapprox(hill_number(c, 1e6), 2.0; rtol = 1e-4)      # → 1 / max p
    sfs = [4, 0, 3]                                   # N = 3, A = 13, k_min = 1
    @test sfs_hill_number(sfs, -Inf) ≈ 13.0           # A / k_min
    @test sfs_hill_number(sfs, -Inf; frequencies = :cell_fraction) ≈ 3.0   # N / k_min
    for f in (:normalised, :cell_fraction)
        @test isapprox(sfs_hill_number(sfs, -1e6; frequencies = f), sfs_hill_number(sfs, -Inf; frequencies = f); rtol = 1e-4)
        @test isapprox(sfs_hill_number(sfs, 1e6; frequencies = f), sfs_hill_number(sfs, Inf; frequencies = f); rtol = 1e-4)
    end
    @test sfs_renyi_entropy([0, 2], -Inf) ≈ log(2 * 2 / 2)       # one class: A / k_min = 4 / 2
end
