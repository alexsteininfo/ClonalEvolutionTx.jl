base = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, response = LinearResponse(4.0))
start(n) = initialize_population(n; mode = :cells)

@testset "extinction mid-schedule" begin   # Review Focus 1 (AtBlockEnd part: Task 12)
    pop = simulate!(start(1000), MTD(base; N_prog = 2000), StableRNG(1))
    @test schedule_outcome(pop) == :extinct && popsize(pop) == 0 && pop.log[end].label == "MTD"
end

@testset "MTD with resistance progresses" begin
    res = remake(base; mutations = (Resistance(ν = 1e-3, α = 0.0),))
    pop = simulate!(start(10_000), MTD(res; N_prog = 12_000), StableRNG(2))
    @test schedule_outcome(pop) in (:progressed, :extinct) && pop.log[end].label == "MTD"
end

@testset "adaptive on/off cycles between thresholds" begin
    pop = start(10_000); t0 = pop.t
    simulate!(pop, AdaptiveOnOff(base; N_low = 5000, N_up = 10_000, N_prog = 20_000, t_max = 30.0), StableRNG(3))
    labels = [e.label for e in pop.log]
    @test labels[1:4] == ["AT on", "AT off", "AT on", "AT off"]
    @test all(e.end_reason == :below for e in pop.log[1:end-1] if e.label == "AT on")
    @test all(e.end_reason == :above for e in pop.log[1:end-1] if e.label == "AT off")
    @test schedule_outcome(pop) == :tmax && pop.t == t0 + 30.0
    @test_throws ArgumentError AdaptiveOnOff(base; N_low = 0, N_up = 10, N_prog = 20, t_max = 1.0)
    @test_throws ArgumentError AdaptiveOnOff(base; N_low = 5, N_up = 10, N_prog = 20, t_max = Inf)
end

@testset "metronomic switches exactly at the interval ends" begin
    mild = remake(base; response = LinearResponse(0.6))     # net −0.1 when on: no extinction
    pop = start(1000); t0 = pop.t
    simulate!(pop, Metronomic(mild; intervals = [(0.0, 1.0), (2.0, 3.0)], N_prog = 10^6, t_max = 4.0), StableRNG(4))
    @test [e.label for e in pop.log] == ["metronomic on", "metronomic off", "metronomic on", "metronomic off"]
    @test [e.t_end for e in pop.log] ≈ t0 .+ [1.0, 2.0, 3.0, 4.0]
end

@testset "stabilising: reduce then hold" begin
    pop = simulate!(start(10_000), Stabilising(base; N_stable = 5000, N_prog = 20_000, t_max = 5.0), StableRNG(5))
    @test pop.log[1].label == "stabilising: reduce" && pop.log[1].end_reason == :below
    @test pop.log[2].label == "stabilising: hold" && pop.log[2].dose == "stabilising"
end

@testset "custom and user schedules" begin
    pop = simulate!(start(1000), CustomDose(base, SizeDose(N -> N > 1200 ? 1.0 : 0.0); N_prog = 10^6, t_max = 2.0), StableRNG(6))
    @test pop.log[end].label == "custom dose" && pop.log[end].dose == "m(N)"
    two = FunctionSchedule((pop, i, t0) -> i < 2 ? remake(base; tmax = t0 + 0.5 * (i + 1), label = "b$i") : nothing)
    pop2 = simulate!(start(10), two, StableRNG(7))
    @test [e.label for e in pop2.log] == ["b0", "b1"]
end
