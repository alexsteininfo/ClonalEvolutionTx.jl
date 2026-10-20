growth(; kw...) = BirthDeathBlock(; birth = s -> 1.0, death = s -> 0.0, kw...)

@testset "trajectory, AtTime, AtSize, AtBlockEnd" begin
    spec = MeasurementSpec(trajectory_dt = 0.5, triggers = [AtTime(1.0), AtSize(500), AtBlockEnd()],
                           stats = [popsize, :clones => p -> length(living_clones(p.clones))])
    acc = MeasurementAccumulator(spec)
    pop = simulate!(initialize_population(1), growth(stop_above = 1000), StableRNG(1); accumulator = acc)
    m = finalize_measurements(acc)
    @test [p.t for p in m.trajectory] == collect(0.0:0.5:pop.t)
    s1, s2, s3 = m.snapshots
    @test s1.trigger == AtTime(1.0) && s1.t == 1.0 && s1.stats.popsize == s1.N
    @test s2.trigger == AtSize(500) && s2.N == 500 && s2.stats.clones == 1
    @test s3.trigger == AtBlockEnd() && s3.N == 1000 && s3.label == "block"
    ref = simulate!(initialize_population(1), growth(tmax = 1.0), StableRNG(1))
    @test s1.N == popsize(ref)          # AtTime sees the state exactly at t
end

@testset "AtTime boundaries" begin   # Review Focus 5
    pop = simulate!(initialize_population(1), growth(tmax = 2.0), StableRNG(2))
    @test_throws ArgumentError simulate!(pop, growth(tmax = 3.0), StableRNG(3);
                                         accumulator = MeasurementAccumulator(MeasurementSpec(triggers = [AtTime(1.0)])))
    acc = MeasurementAccumulator(MeasurementSpec(triggers = [AtTime(4.0)], stats = [popsize]))
    simulate!(pop, growth(tmax = 4.0), StableRNG(4); accumulator = acc)
    simulate!(pop, growth(tmax = 5.0), StableRNG(5); accumulator = acc)
    snaps = finalize_measurements(acc).snapshots
    @test count(s -> s.trigger == AtTime(4.0), snaps) == 1 && snaps[1].t == 4.0
end

@testset "AtSize down, every crossing, When" begin
    pop0 = initialize_population(200; mode = :cells)
    shrink = BirthDeathBlock(birth = s -> 0.0, death = s -> 1.0, stop_below = 50)
    acc = MeasurementAccumulator(MeasurementSpec(triggers = [AtSize(100; direction = :down), When(p -> popsize(p) == 150)]))
    simulate!(pop0, shrink, StableRNG(6); accumulator = acc)
    snaps = finalize_measurements(acc).snapshots
    @test [s.N for s in snaps] == [150, 100]
    critical = BirthDeathBlock(birth = s -> 1.0, death = s -> 1.0, tmax = 10.0)
    acc2 = MeasurementAccumulator(MeasurementSpec(triggers = [AtSize(120; direction = :either, every = true)]))
    simulate!(initialize_population(120; mode = :cells), critical, StableRNG(12); accumulator = acc2)
    snaps2 = finalize_measurements(acc2).snapshots
    @test length(snaps2) > 1 && all(s -> s.N == 120, snaps2)
end

@testset "trajectory clone sizes are opt-in" begin
    two = [(n = 10, s = 0.0, α = 1.0), (n = 5, s = 0.0, α = 0.0)]
    acc = MeasurementAccumulator(MeasurementSpec(trajectory_dt = 1.0, trajectory_clone_sizes = true))
    simulate!(initialize_population(two), growth(tmax = 2.0), StableRNG(13); accumulator = acc)
    tr = finalize_measurements(acc).trajectory
    @test tr[1].clone_sizes == [10, 5] && tr[1].n_sensitive == 10 && tr[1].n_resistant == 5
    acc0 = MeasurementAccumulator(MeasurementSpec(trajectory_dt = 1.0))
    simulate!(initialize_population(two), growth(tmax = 2.0), StableRNG(13); accumulator = acc0)
    @test finalize_measurements(acc0).trajectory[1].clone_sizes === nothing
end

@testset "AtBlockEnd fires when a schedule ends in extinction" begin   # Review Focus 1
    base = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, response = LinearResponse(4.0))
    acc = MeasurementAccumulator(MeasurementSpec(triggers = [AtBlockEnd()], stats = [popsize]))
    pop = simulate!(initialize_population(1000; mode = :cells), MTD(base; N_prog = 2000), StableRNG(14); accumulator = acc)
    @test schedule_outcome(pop) == :extinct && finalize_measurements(acc).snapshots[end].N == 0
end

@testset "samples at extinction are recorded as nothing" begin   # final review F1
    base = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, response = LinearResponse(4.0))
    acc = MeasurementAccumulator(MeasurementSpec(triggers = [AtBlockEnd()], stats = [popsize],
                                                 samples = [50, 100]))
    pop = simulate!(initialize_population(1000; mode = :cells), MTD(base; N_prog = 2000), StableRNG(14); accumulator = acc)
    sn = finalize_measurements(acc).snapshots[end]
    @test sn.N == 0 && sn.samples == [(n = 50, stats = nothing), (n = 100, stats = nothing)]
    @test pop.log[end].end_reason == :extinct && schedule_outcome(pop) == :extinct
    @test sn.block == length(pop.log)
end

@testset "a sample larger than the population is nothing, the others are drawn" begin   # F1
    pop = initialize_population(75)
    acc = MeasurementAccumulator(MeasurementSpec(triggers = [AtBlockEnd()], stats = [popsize, cell_ids],
                                                 samples = [50, 100, 75], sample_seed = 3))
    simulate!(pop, BirthDeathBlock(birth = s -> 0.0, death = s -> 0.0, tmax = 1.0), StableRNG(1); accumulator = acc)
    s = only(finalize_measurements(acc).snapshots).samples
    @test [x.n for x in s] == [50, 100, 75] && s[2] == (n = 100, stats = nothing)
    for k in (1, 3)                     # seeds as without the skipped size: derive_seed(seed, index, n)
        n = s[k].n
        @test s[k].stats == measure(sample_population(pop, n; seed = derive_seed(3, 1, n)), [popsize, cell_ids])
    end
end

@testset "measure never changes the population" begin
    pop = simulate!(initialize_population(1), growth(stop_above = 100), StableRNG(7))
    ids = copy(cell_ids(pop))
    r = measure(pop, [popsize, :burden => p -> sum(mutations_per_cell(p))])
    @test r.popsize == 100 && cell_ids(pop) == ids
end

@testset "restarts roll the accumulator back" begin
    block = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.9, stop_above = 50, restart_on_extinction = true)
    acc = MeasurementAccumulator(MeasurementSpec(trajectory_dt = 0.1))
    pop = simulate!(initialize_population(1), block, StableRNG(11); accumulator = acc)
    ts = [p.t for p in finalize_measurements(acc).trajectory]
    @test issorted(ts) && allunique(ts) && ts[end] <= pop.t
end

# ── Beyond the brief: the controller's rulings 6, 8 and 9, review fixes ───────

@testset "statistics see the clock at the measurement time" begin   # SPEC §7.1
    clock(p) = p.t
    spec = MeasurementSpec(trajectory_dt = 0.5, triggers = [AtTime(1.0), Every(0.75), AtTime(2.0)],
                           stats = [clock, popsize])
    acc = MeasurementAccumulator(spec)
    pop = simulate!(initialize_population(1), growth(tmax = 2.0), StableRNG(1); accumulator = acc)
    snaps = finalize_measurements(acc).snapshots
    @test [s.t for s in snaps] == [0.0, 0.75, 1.0, 1.5, 2.0]          # AtTime(2.0) at the :tmax end
    @test all(s -> s.stats.clock == s.t, snaps)
    ref = simulate!(initialize_population(1), growth(tmax = 2.0), StableRNG(1))   # the run is unchanged
    @test pop.t == ref.t && pop.log == ref.log && cell_ids(pop) == cell_ids(ref) && pop.tree.birthtime == ref.tree.birthtime
    @test snaps[3].N == popsize(simulate!(initialize_population(1), growth(tmax = 1.0), StableRNG(1)))
end

@testset "firing does not repeat across chained blocks" begin   # Review Focus 5, Ruling 6
    acc = MeasurementAccumulator(MeasurementSpec(trajectory_dt = 0.5,
        triggers = [AtSize(50; every = true), When(p -> popsize(p) == 50; every = true), Every(0.5)]))
    pop = simulate!(initialize_population(1), growth(stop_above = 50), StableRNG(8); accumulator = acc)
    simulate!(pop, growth(tmax = pop.t), StableRNG(9); accumulator = acc)          # ends at once
    simulate!(pop, growth(stop_above = 60), StableRNG(10); accumulator = acc)
    m = finalize_measurements(acc)
    @test [s.block for s in m.snapshots if s.trigger isa Union{AtSize,When}] == [1, 1]
    grid = collect(0.0:0.5:pop.t)
    @test [p.t for p in m.trajectory] == grid
    @test [s.t for s in m.snapshots if s.trigger == Every(0.5)] == grid
end

@testset "AtSize needs an arrival, When fires at block start" begin   # Ruling 6
    acc = MeasurementAccumulator(MeasurementSpec(triggers = [AtSize(10), When(p -> popsize(p) >= 10)]))
    simulate!(initialize_population(10), growth(stop_above = 12), StableRNG(15); accumulator = acc)
    snaps = finalize_measurements(acc).snapshots
    @test length(snaps) == 1 && snaps[1].trigger isa When && snaps[1].t == 0.0
end

@testset "the recorded level is the dose at the measurement time" begin   # Ruling 8
    pop = simulate!(initialize_population(100; mode = :cells), growth(tmax = 1.0), StableRNG(16))
    ramp = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.2, dose = TimeDose(t -> t / 4), tmax = 3.0)
    acc = MeasurementAccumulator(MeasurementSpec(trajectory_dt = 0.5, triggers = [AtTime(2.0), AtBlockEnd()]))
    simulate!(pop, ramp, StableRNG(17); accumulator = acc)
    m = finalize_measurements(acc)
    @test [p.t for p in m.trajectory] == 1.0:0.5:3.0                 # the grid starts at the accumulator's start
    @test [p.m for p in m.trajectory] == [0.0, 0.125, 0.25, 0.375, 0.5]
    @test [s.m for s in m.snapshots] == [0.25, 0.5]
end

@testset "a time grid skips what the accumulator did not see" begin   # Ruling 9
    acc = MeasurementAccumulator(MeasurementSpec(trajectory_dt = 1.0, triggers = [AtTime(5.0)]))
    pop = simulate!(initialize_population(5), growth(tmax = 1.5), StableRNG(18); accumulator = acc)
    simulate!(pop, growth(tmax = 3.5), StableRNG(19))                      # not measured
    simulate!(pop, growth(tmax = 5.0), StableRNG(20); accumulator = acc)
    m = finalize_measurements(acc)
    @test [p.t for p in m.trajectory] == [0.0, 1.0, 4.0, 5.0]
    @test length(m.snapshots) == 1 && m.snapshots[1].t == 5.0 && m.snapshots[1].block == 3
end

@testset "spec and statistic validation" begin
    @test_throws ArgumentError MeasurementSpec(trajectory_dt = 0.0)
    @test_throws ArgumentError MeasurementSpec(trajectory_clone_sizes = true)
    @test_throws ArgumentError MeasurementSpec(stats = [popsize, popsize])
    @test_throws ArgumentError MeasurementSpec(stats = [p -> 1])               # anonymous: name it
    @test_throws ArgumentError MeasurementSpec(stats = [42])
    @test_throws ArgumentError MeasurementSpec(samples = [0])
    @test_throws ArgumentError AtSize(10; direction = :sideways)
    @test_throws ArgumentError AtTime(NaN)
    @test_throws ArgumentError Every(0.0)
    @test AtSize(10) == AtSize(10; direction = :up, every = false) && AtSize(10) != AtSize(10; every = true)
    @test measure(initialize_population(3), []) == NamedTuple()
end
