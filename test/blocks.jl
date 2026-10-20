@testset "block construction" begin
    b = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5)
    @test b.μ == 0.0 && b.tmax == Inf && b.mutations == () && b.label == "block"
    b2 = remake(b; μ = 2.0, label = "growth")
    @test b2.μ == 2.0 && b2.label == "growth" && b2.death(0.0) == 0.5
    @test_throws ArgumentError BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, μ = -1.0)
    @test_throws ArgumentError BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, tmax = NaN)
    @test_throws ArgumentError BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, max_clones = 0)
    @test_throws ArgumentError BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, mutations = (1.0,))
    @test_throws ArgumentError BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, mutations = [Drivers(ν = 0.1, Δs = 0.1), "x"])
    v = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, mutations = [Drivers(ν = 0.1, Δs = 0.1)])
    @test v.mutations isa Tuple && length(v.mutations) == 1
    @test_throws ArgumentError remake(b; nonsense = 1)
end

@testset "divide and die keep both modes in step" begin
    block = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.0, μ = 2.0)
    pops = map((:tree, :cells)) do mode
        pop = initialize_population(3; mode); rng = StableRNG(1)
        ClonalEvolutionTx._start_block!(pop, block, rng)
        ClonalEvolutionTx._divide!(pop, 2, 1.0, block, rng)
        ClonalEvolutionTx._die!(pop, 1)
        pop
    end
    @test [popsize(p) for p in pops] == [3, 3]
    @test cell_ids(pops[1]) == cell_ids(pops[2])
    @test mutations_per_cell(pops[1]) == mutations_per_cell(pops[2])
    @test divisions_per_cell(pops[1]) == divisions_per_cell(pops[2])
end

@testset "mutations nest in list order" begin
    block = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.0,
                            mutations = (Drivers(ν = 1.0, Δs = 0.1), Resistance(ν = 1.0, α = 0.5)))
    pop = initialize_population(1); rng = StableRNG(2)
    ClonalEvolutionTx._start_block!(pop, block, rng)
    ClonalEvolutionTx._divide!(pop, 1, 1.0, block, rng)
    @test n_clones(pop.clones) == 5
    @test pop.clones.kind[2:5] == [:driver, :resistance, :driver, :resistance]
    @test pop.clones.parent[2:5] == [1, 2, 1, 4]
    @test living_clones(pop.clones) == [3, 5] && pop.clones.size[1] == 0
    @test traits(pop.clones, 3) == Traits(0.1, 0.5)
end

@testset "max_clones caps new clones" begin
    block = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.0,
                            mutations = (Drivers(ν = 1.0, Δs = 0.1),), max_clones = 2)
    pop = initialize_population(1); rng = StableRNG(3)
    ClonalEvolutionTx._start_block!(pop, block, rng)
    ClonalEvolutionTx._divide!(pop, 1, 1.0, block, rng)
    @test n_clones(pop.clones) == 2
end

@testset "neutral draws never touch the dynamics stream" begin
    states = map((0.0, 3.0)) do μ
        block = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.0, μ = μ)
        pop = initialize_population(1); rng = StableRNG(4)
        ClonalEvolutionTx._start_block!(pop, block, rng)
        ClonalEvolutionTx._divide!(pop, 1, 1.0, block, rng)
        rand(rng)
    end
    @test states[1] == states[2]
end

@testset "size-1 clone dividing without mutation keeps its clone" begin
    block = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.0)
    for mode in (:tree, :cells)
        pop = initialize_population(1; mode); rng = StableRNG(5)
        ClonalEvolutionTx._start_block!(pop, block, rng)
        ClonalEvolutionTx._divide!(pop, 1, 1.0, block, rng)
        @test pop.clones.size[1] == 2 && popsize(pop) == 2
    end
end

@testset "tree compaction in _die! keeps living cells" begin
    block = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.0, μ = 1.0)
    pops = map((:tree, :cells)) do mode
        pop = initialize_population(2; mode); rng = StableRNG(6)
        ClonalEvolutionTx._start_block!(pop, block, rng)
        for i in 1:4200
            ClonalEvolutionTx._divide!(pop, 1 + (7i) % popsize(pop), Float64(i), block, rng)
        end
        mode === :tree && (@test length(pop.tree) > 8000)
        while popsize(pop) > 5     # mass die-off: dead rows pass the compaction trigger
            ClonalEvolutionTx._die!(pop, popsize(pop))
        end
        pop
    end
    @test length(pops[1].tree) < 1000      # compaction ran
    @test cell_ids(pops[1]) == cell_ids(pops[2])
    @test mutations_per_cell(pops[1]) == mutations_per_cell(pops[2])
    @test divisions_per_cell(pops[1]) == divisions_per_cell(pops[2])
    @test all(i -> ClonalEvolutionTx._isalive(pops[1].tree, i), pops[1].alive)
    @test sum(pops[1].clones.size) == popsize(pops[1])
end

@testset "divide and die do not allocate" begin
    block = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.0, μ = 2.0,
                            mutations = (Drivers(ν = 0.0, Δs = 0.1), Resistance(ν = 0.0, α = 0.5)))
    for mode in (:tree, :cells)
        pop = initialize_population(2; mode); rng = StableRNG(7)
        ClonalEvolutionTx._start_block!(pop, block, rng)
        for i in 1:200   # grow arrays past reallocation
            ClonalEvolutionTx._divide!(pop, 1, 1.0, block, rng)
        end
        step!(pop) = (ClonalEvolutionTx._divide!(pop, 1, 1.0, block, rng); ClonalEvolutionTx._die!(pop, 2); nothing)
        step!(pop)
        @test (@allocated step!(pop)) == 0
    end
end
