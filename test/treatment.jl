v1(; kw...) = BirthDeathBlock(; birth = s -> 1.0, death = s -> 0.5, response = LinearResponse(4.0), kw...)
meanN(pop0, block; runs, seed0 = 0) = mean(popsize(simulate!(copy(pop0), block, StableRNG(seed0 + s))) for s in 1:runs)

@testset "a fully resistant clone ignores the dose" begin
    p0 = initialize_population(1000; mode = :cells, α = 0.0)
    @test isapprox(meanN(p0, v1(dose = FixedDose(1.0), tmax = 2.0); runs = 200), 1000 * exp(1.0); rtol = 0.03)
end

@testset "full dose kills sensitive cells at rate b - d0 - d_S" begin
    p0 = initialize_population(10_000; mode = :cells)
    @test isapprox(meanN(p0, v1(dose = FixedDose(1.0), tmax = 0.5); runs = 100), 10_000 * exp(-1.75); rtol = 0.03)
end

@testset "stabilising dose holds the population" begin
    p0 = initialize_population(2000; mode = :cells)
    @test isapprox(meanN(p0, v1(dose = StabilisingDose(), tmax = 20.0); runs = 100), 2000; rtol = 0.05)
end

@testset "stabilising dose with only fully resistant cells" begin
    pop = simulate!(initialize_population(100; mode = :cells, α = 0.0), v1(dose = StabilisingDose(), tmax = 1.0), StableRNG(1))
    @test pop.m == 1.0 && pop.log[end].dose == "stabilising"
end

@testset "time dose: constant equals fixed, step equals the analytic mean" begin
    p0 = initialize_population(1000; mode = :cells)
    @test isapprox(meanN(p0, v1(dose = TimeDose(t -> 0.5), tmax = 1.0); runs = 300), 1000 * exp(-1.5); rtol = 0.03)
    step = TimeDose(t -> t < 0.5 ? 1.0 : 0.0)
    @test isapprox(meanN(p0, v1(dose = step, tmax = 1.0); runs = 300, seed0 = 500), 1000 * exp(-1.5); rtol = 0.03)
    @test_throws ArgumentError simulate!(copy(p0), v1(dose = TimeDose(t -> 0.9; m_max = 0.5), tmax = 1.0), StableRNG(2))
end

@testset "size dose follows N" begin
    pop = simulate!(initialize_population(1000; mode = :cells), v1(dose = SizeDose(N -> N >= 1500 ? 1.0 : 0.0), tmax = 10.0), StableRNG(3))
    @test popsize(pop) < 2000 && pop.m == (popsize(pop) >= 1500 ? 1.0 : 0.0)
end

@testset "stabilising dose once only resistant cells remain (running-sum round-off)" begin
    p0 = initialize_population([(n = 200, s = 0.0, α = 0.1), (n = 200, s = 0.0, α = 0.7), (n = 20, s = 0.0, α = 0.0)]; mode = :cells)
    for seed in 1:10
        pop = simulate!(copy(p0), v1(response = LinearResponse(20.0), dose = StabilisingDose(), tmax = 30.0, stop_above = 50_000), StableRNG(seed))
        @test pop.clones.maxα == 0.0          # precondition: no α > 0 clone alive
        @test pop.m == 1.0 && pop.log[end].m_end == 1.0
    end
end

# ── Final review fixes ────────────────────────────────────────────────────────

# Not monotone: D(0.5) = 1 exceeds D(1) = 0.5.
struct BumpResponse <: AbstractResponse end
ClonalEvolutionTx.drug_death(::BumpResponse, m::Real) = m <= 0.5 ? 2m : 1.5 - m

@testset "time dose with a Hill response near m_max" begin   # final review F2
    hill = HillResponse(1.0, 3.5, 0.5)
    @test drug_death(hill, prevfloat(1.0, 2)) > drug_death(hill, 1.0)   # precondition: rounding
    # a monotone response rounded a few ulps above D(m_max) is no error
    near = v1(response = hill, dose = TimeDose(t -> prevfloat(1.0, 2)), tmax = 1.0)
    @test popsize(simulate!(initialize_population(100; mode = :cells), near, StableRNG(1))) > 0
    # the review's reproduction: a dose approaching m_max, run to extinction
    approach = v1(response = hill, dose = TimeDose(t -> 1 - exp(-(t + 20))))
    pops = [simulate!(initialize_population(20_000; mode = :cells), approach, StableRNG(s)) for s in 1:5]
    @test all(p -> p.log[end].end_reason == :extinct, pops)
    # a response that is not monotone is still an error
    bump = v1(response = BumpResponse(), dose = TimeDose(t -> 0.5), tmax = 1.0)
    @test_throws ArgumentError simulate!(initialize_population(100; mode = :cells), bump, StableRNG(2))
    # with D(m_max) = 0 and no other rate no event can happen: still detected, never a hang
    idle = BirthDeathBlock(birth = s -> 0.0, death = s -> 0.0, response = LinearResponse(0.0),
                           dose = TimeDose(t -> 0.5))
    @test_throws ErrorException simulate!(initialize_population(4), idle, StableRNG(3))
    @test simulate!(initialize_population(4), remake(idle; tmax = 2.0), StableRNG(3)).t == 2.0
end

@testset "a time dose ends at its level at t_end" begin   # final review m1
    pop = simulate!(initialize_population(100; mode = :cells), v1(tmax = 1.0), StableRNG(4))
    t_start = pop.t
    simulate!(pop, v1(dose = TimeDose(t -> t / 4), tmax = 3.0), StableRNG(5))
    @test pop.log[end].end_reason == :tmax
    @test pop.log[end].m_end == (3.0 - t_start) / 4 && pop.m == (3.0 - t_start) / 4
end
