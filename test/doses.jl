@testset "responses and inverses" begin
    @test drug_death(LinearResponse(4.0), 0.5) == 2.0
    h = HillResponse(4.0, 2.0, 0.5)
    @test drug_death(h, 0.5) ≈ 2.0
    for r in (LinearResponse(4.0), h), m in (0.0, 0.1, 0.5, 0.9, 1.0)
        @test inverse_drug_death(r, drug_death(r, m)) ≈ m atol = 1e-10
    end
    @test_throws ArgumentError LinearResponse(-1.0)
    @test_throws ArgumentError HillResponse(0.0, 2.0, 0.5)
    @test_throws ArgumentError HillResponse(4.0, -2.0, 0.5)
    @test_throws ArgumentError HillResponse(4.0, 2.0, 0.0)
    @test_throws ArgumentError inverse_drug_death(LinearResponse(4.0), 5.0)
    @test_throws ArgumentError inverse_drug_death(LinearResponse(4.0), -0.1)
    @test_throws ArgumentError inverse_drug_death(h, 3.9)
    @test inverse_drug_death(LinearResponse(0.0), 0.0) == 0.0
    @test_throws ArgumentError inverse_drug_death(LinearResponse(0.0), 0.1)
end

struct SquareResponse <: AbstractResponse end
ClonalEvolutionTx.drug_death(::SquareResponse, m) = 3m^2
@testset "custom response falls back to bisection" begin
    @test inverse_drug_death(SquareResponse(), 0.75) ≈ 0.5 atol = 1e-9
    @test_throws ArgumentError inverse_drug_death(SquareResponse(), 3.5)
end

@testset "stabilising level reproduces manuscript_v1 eq. 7" begin
    b, d0, l, N, S, R, α = 1.0, 0.5, 1e-4, 1000, 700, 300, 0.2
    m_v1 = (b - d0 - l * N) * N / (4.0 * (S + α * R))
    @test ClonalEvolutionTx._stabilising_level(LinearResponse(4.0), 1.0; B = 1.0, C = l * N, N = N,
        sumNβ = b * N, sumNδ = d0 * N, sumNα = S + α * R) ≈ m_v1
    h = HillResponse(4.0, 2.0, 0.5); x = 1.5          # target D(m) = 1.5
    @test ClonalEvolutionTx._stabilising_level(h, 1.0; B = 1.0, C = 0.0, N = 100,
        sumNβ = 100.0 + 150.0, sumNδ = 100.0, sumNα = 100.0) ≈ 0.5 * (x / (4.0 - x))^(1 / 2)
end

@testset "stabilising level edge cases" begin   # Review Focus 4
    lvl(; sumNδ = 50.0, sumNα = 100.0, m_max = 1.0) = ClonalEvolutionTx._stabilising_level(
        LinearResponse(4.0), m_max; B = 1.0, C = 0.0, N = 100, sumNβ = 100.0, sumNδ, sumNα)
    @test lvl(sumNδ = 200.0) == 0.0          # population shrinks without drug
    @test lvl(sumNα = 0.0) == 1.0            # no sensitive cells: m = m_max, no NaN
    @test lvl(sumNα = -1e-13) == 1.0         # round-off residue of a running sum
    @test lvl(sumNα = 1.0) == 1.0            # target above D(m_max)
    @test lvl(sumNα = 1.0, m_max = 0.6) == 0.6
end

@testset "dose validation and levels" begin
    @test_throws ArgumentError FixedDose(1.5)
    @test_throws ArgumentError FixedDose(-0.1)
    @test_throws ArgumentError StabilisingDose(m_max = 0.0)
    @test_throws ArgumentError StabilisingDose(m_max = 1.5)
    @test_throws ArgumentError TimeDose(t -> 0.1; m_max = 0.0)
    @test ClonalEvolutionTx._size_level(SizeDose(N -> N > 100 ? 1.0 : 0.0), 150) == 1.0
    @test_throws ArgumentError ClonalEvolutionTx._size_level(SizeDose(N -> 2.0), 10)
    @test ClonalEvolutionTx._time_level(TimeDose(t -> 0.4; m_max = 0.5), 3.0) == 0.4
    @test_throws ArgumentError ClonalEvolutionTx._time_level(TimeDose(t -> 0.8; m_max = 0.5), 0.0)
    @test_throws ArgumentError ClonalEvolutionTx._time_level(TimeDose(t -> -0.1), 0.0)
    @test ClonalEvolutionTx._describe(FixedDose(1.0)) == "fixed m = 1.0"
    @test ClonalEvolutionTx._describe(StabilisingDose()) == "stabilising"
    @test ClonalEvolutionTx._describe(SizeDose(identity)) == "m(N)"
    @test ClonalEvolutionTx._describe(TimeDose(identity)) == "m(t)"
end

@testset "derive_seed" begin
    @test derive_seed(1, 2) == derive_seed(1, 2)
    @test derive_seed(1, 2) != derive_seed(1, 3) != derive_seed(2, 2)
    @test derive_seed(1, 2) isa UInt64
    @test derive_seed(1) != derive_seed(1, 0)
end
