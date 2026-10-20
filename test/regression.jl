# Golden-value regression tests.
#
# These freeze the outputs of two fixed-seed runs (StableRNG), so that any change in how
# the simulator consumes random numbers is caught: stored simulation results must stay
# reproducible from their seeds. Any change in random-number consumption changes these
# values and needs the author's sign-off; regenerate them only deliberately.

@testset "golden: neutral growth, tree mode (StableRNG(2026))" begin
    rng = StableRNG(2026)
    pop = initialize_population(1)
    block = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, μ = 1.0,
                            stop_above = 5000, restart_on_extinction = true)
    simulate!(pop, block, rng)

    @test popsize(pop) == 5000
    @test pop.t ≈ 17.775402387665398 rtol = 1e-10
    @test n_clones(pop.clones) == 1
    @test sum(mutations_per_cell(pop)) == 166530
    @test site_frequency_spectrum(pop)[1:5] == [6146, 2242, 1263, 826, 542]
    @test clone_sizes(pop)[1:min(end, 5)] == [5000]
    @test length(pop.log) == 1
end

@testset "golden: resistance growth then AdaptiveOnOff, cell mode (StableRNG(2027))" begin
    rng = StableRNG(2027)
    base = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5,
                           death_density = N -> 0.5 / 1e5 * N,
                           response = LinearResponse(4.0), μ = 1.0,
                           mutations = (Resistance(ν = 0.002, α = 0.8),))
    pop = initialize_population(1; mode = :cells)
    simulate!(pop, remake(base; stop_above = 5000, restart_on_extinction = true), rng)
    @test popsize(pop) == 5000
    @test pop.t ≈ 14.185377275085791 rtol = 1e-10

    simulate!(pop, AdaptiveOnOff(base; N_low = 2500, N_up = 5000, N_prog = 10_000,
                                 t_max = 20.0), rng)
    @test popsize(pop) == 4717
    @test pop.t ≈ 34.18537727508579 rtol = 1e-10
    @test n_clones(pop.clones) == 318
    @test sum(mutations_per_cell(pop)) == 320249
    @test clone_sizes(pop)[1:min(end, 5)] == [2086, 408, 1178, 77, 7]
    @test length(pop.log) == 25
    @test schedule_outcome(pop) == :tmax
end
