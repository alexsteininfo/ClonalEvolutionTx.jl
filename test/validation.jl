@testset "extinction probability d/b" begin
    block = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, stop_above = 30)
    ext = mean(popsize(simulate!(initialize_population(1), block, StableRNG(s))) == 0 for s in 1:2000)
    @test abs(ext - 0.5) < 0.04
end

@testset "E N(t) = N0 exp((b - d) t)" begin
    block = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, tmax = 2.0)
    m = mean(popsize(simulate!(initialize_population(100; mode = :cells), block, StableRNG(s))) for s in 1:400)
    @test isapprox(m, 100 * exp(1.0); rtol = 0.03)
end

@testset "logistic growth settles at K" begin
    block = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.5, death_density = N -> 0.5 / 1000 * N, tmax = 20.0)
    m = mean(popsize(simulate!(initialize_population(1000; mode = :cells), block, StableRNG(s))) for s in 1:100)
    @test isapprox(m, 1000; rtol = 0.05)
end

@testset "mean burden = μ × mean divisions" begin
    ratios = map(1:40) do s
        pop = grown(:cells; N = 10_000, seed = s, μ = 2.0)
        mean(mutations_per_cell(pop)) / mean(divisions_per_cell(pop))
    end
    @test isapprox(mean(ratios), 2.0; rtol = 0.03)
end

@testset "resistant clones created ≈ 2ν × divisions" begin
    block = BirthDeathBlock(birth = s -> 1.0, death = s -> 0.0, stop_above = 10_000,
                            mutations = (Resistance(ν = 0.01, α = 0.5),))
    created = mean(n_clones(simulate!(initialize_population(1; mode = :cells), block, StableRNG(s)).clones) - 1 for s in 1:20)
    @test isapprox(created, 2 * 0.01 * 9_999; rtol = 0.05)
end
