pure_birth(; kw...) = BirthDeathBlock(; birth = s -> 1.0, death = s -> 0.0, kw...)

@testset "pure birth reaches the stop size" for mode in (:tree, :cells)
    pop = simulate!(initialize_population(1; mode), pure_birth(stop_above = 1000), StableRNG(1))
    @test popsize(pop) == 1000 && pop.log[end].end_reason == :above
    @test pop.log[end].ticks == pop.log[end].accepted == 999      # uniform rates: nothing rejected
end

@testset "tmax is exact" begin
    pop = simulate!(initialize_population(1), pure_birth(tmax = 2.5), StableRNG(2))
    @test pop.t == 2.5 && pop.log[end].end_reason == :tmax
    @test all(<=(2.5), pop.tree.birthtime)
end

@testset "extinction and restarts" begin
    pop = simulate!(initialize_population(10), BirthDeathBlock(birth = s -> 0.0, death = s -> 1.0), StableRNG(3))
    @test popsize(pop) == 0 && pop.log[end].end_reason == :extinct && isfinite(pop.t)
    block = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.9, stop_above = 100, restart_on_extinction = true)
    runs = [simulate!(initialize_population(1), block, StableRNG(s)) for s in 1:20]
    @test all(p -> popsize(p) == 100, runs) && any(p -> p.log[end].restarts > 0, runs)
end

@testset "zero total rate" begin   # Review Focus 2
    idle = BirthDeathBlock(birth = s -> 0.0, death = s -> 0.0, tmax = 5.0)
    pop = simulate!(initialize_population(4), idle, StableRNG(4))
    @test pop.t == 5.0 && popsize(pop) == 4 && pop.log[end].end_reason == :tmax
    @test_throws ErrorException simulate!(initialize_population(4), remake(idle; tmax = Inf), StableRNG(4))
end

@testset "stop conditions are checked at block start" begin
    pop = initialize_population(50)
    simulate!(pop, pure_birth(stop_above = 10), StableRNG(5))
    @test popsize(pop) == 50 && pop.log[end].ticks == 0 && pop.log[end].end_reason == :above
end

@testset "chaining continues clock, tree and log" begin
    pop = initialize_population(1)
    simulate!(pop, pure_birth(stop_above = 100), StableRNG(6)); t1 = pop.t
    simulate!(pop, remake(pure_birth(); tmax = t1 + 1.0, label = "second"), StableRNG(7))
    @test pop.t == t1 + 1.0 && length(pop.log) == 2 && pop.log[2].label == "second" && pop.log[2].t_start == t1
end

@testset "extinct clones leave the bound" begin
    # one sensitive cell among 1000 fully resistant ones under full dose
    pop = initialize_population([(n = 1, s = 0.0, α = 1.0), (n = 1000, s = 0.0, α = 0.0)])
    treat = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, response = LinearResponse(4.0),
                            dose = FixedDose(1.0), stop = p -> p.clones.size[1] == 0)
    simulate!(pop, treat, StableRNG(8))
    simulate!(pop, remake(treat; stop = Returns(false), tmax = pop.t + 0.5), StableRNG(9))
    @test pop.log[end].accepted / pop.log[end].ticks > 0.99
end

@testset "log describes the dose" begin
    pop = simulate!(initialize_population(5), remake(pure_birth(); tmax = 0.1, dose = FixedDose(1.0)), StableRNG(10))
    @test pop.log[end].dose == "fixed m = 1.0" && pop.log[end].m_end == 1.0
end

@testset "negative rates are errors, not clamps" begin
    @test_throws ArgumentError simulate!(initialize_population(5),
        BirthDeathBlock(birth = s -> 1.0, death = s -> 0.0, death_density = N -> -1.0, tmax = 1.0), StableRNG(11))
    @test_throws ArgumentError simulate!(initialize_population(5),
        BirthDeathBlock(birth = s -> -1.0, death = s -> 0.0, tmax = 1.0), StableRNG(12))
end

@testset "a restart from an empty population is an error" begin
    pop = simulate!(initialize_population(3), BirthDeathBlock(birth = s -> 0.0, death = s -> 1.0), StableRNG(13))
    @test popsize(pop) == 0
    @test_throws ArgumentError simulate!(pop, remake(pure_birth(); restart_on_extinction = true), StableRNG(14))
    simulate!(pop, pure_birth(), StableRNG(15))     # without restarts: ends at once
    @test pop.log[end].end_reason == :extinct && pop.log[end].ticks == 0
end

@testset "tree and cell mode agree" begin
    block = BirthDeathBlock(birth = s -> 1.0 + s, death = s -> 0.4, μ = 1.0,
                            mutations = (Drivers(ν = 0.05, Δs = 0.2),), stop_above = 300,
                            restart_on_extinction = true)
    pops = [simulate!(initialize_population(1; mode), block, StableRNG(16)) for mode in (:tree, :cells)]
    @test cell_ids(pops[1]) == cell_ids(pops[2]) && cell_clones(pops[1]) == cell_clones(pops[2])
    @test mutations_per_cell(pops[1]) == mutations_per_cell(pops[2])
    @test pops[1].t == pops[2].t && pops[1].log == pops[2].log
end

@testset "grown fixture" for mode in (:tree, :cells)
    pop = grown(mode; N = 200)
    @test popsize(pop) == 200 && pop.log[end].end_reason == :above
    @test sum(pop.clones.size) == 200
end

@testset "tmax before the population's time is an error" begin
    pop = simulate!(initialize_population(1), pure_birth(tmax = 1.0), StableRNG(17))
    @test_throws ArgumentError simulate!(pop, pure_birth(tmax = 0.5), StableRNG(18))
    simulate!(pop, pure_birth(tmax = 1.0), StableRNG(19))   # tmax == t: ends at once
    @test pop.t == 1.0 && pop.log[end].end_reason == :tmax && pop.log[end].ticks == 0
end
