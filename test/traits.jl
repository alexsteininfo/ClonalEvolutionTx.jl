@testset "Traits validate α" begin
    @test Traits(0.1, 0.5).α == 0.5
    @test_throws ArgumentError Traits(0.0, -0.1)
    @test_throws ArgumentError Traits(0.0, 1.2)
end

@testset "Resistance: min rule, cost of either sign" begin
    rng = StableRNG(1)
    parent = Traits(0.0, 0.5)
    child = mutate(Resistance(ν = 1.0, α = 0.8, Δs = -0.1), parent, rng)
    @test child.α == 0.5            # min(0.8, 0.5): never less resistant
    @test child.s ≈ -0.1            # cost applied although α did not change
    @test mutate(Resistance(ν = 1.0, α = 0.2, Δs = 0.05), parent, rng) == Traits(0.05, 0.2)
    @test mutate(Resistance(ν = 1.0, α = 0.3, Δs = α -> -0.1 * (1 - α)), Traits(0.0, 1.0), rng).s ≈ -0.07
end

@testset "Drivers: update rules, α untouched" begin
    rng = StableRNG(2)
    @test mutate(Drivers(ν = 1.0, Δs = 0.1), Traits(0.2, 0.7), rng) == Traits(0.2 + 0.1, 0.7)
    @test mutate(Drivers(ν = 1.0, Δs = 0.1, update = :multiplicative), Traits(0.2, 1.0), rng).s ≈ 1.2 * 1.1 - 1
    @test mutate(Drivers(ν = 1.0, Δs = 0.1, update = (s, Δ) -> max(s, Δ)), Traits(0.2, 1.0), rng).s == 0.2
end

@testset "distributions are drawn per mutation" begin
    rng = StableRNG(3)
    αs = [mutate(Resistance(ν = 1.0, α = Uniform(0, 1)), Traits(0.0, 1.0), rng).α for _ in 1:2000]
    @test 0.47 < mean(αs) < 0.53
    @test_throws ArgumentError [mutate(Resistance(ν = 1.0, α = Normal(0.5, 2.0)), Traits(0.0, 1.0), rng) for _ in 1:100]
end

@testset "constructor validation and kinds" begin
    @test_throws ArgumentError Drivers(ν = -0.1, Δs = 0.1)
    @test_throws ArgumentError Drivers(ν = 1.5, Δs = 0.1)
    @test_throws ArgumentError Resistance(ν = 0.1, α = 1.2)
    @test_throws ArgumentError Drivers(ν = 0.1, Δs = 0.1, update = :foo)
    c = CustomMutation(ν = 0.5, f = (rng, t) -> Traits(t.s + 1, t.α / 2))
    @test mutate(c, Traits(0.0, 1.0), StableRNG(1)) == Traits(1.0, 0.5)
    @test mutation_kind.((Drivers(ν = 0.1, Δs = 0.1), Resistance(ν = 0.1, α = 0.5), c)) == (:driver, :resistance, :custom)
end
