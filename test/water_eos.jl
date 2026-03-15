# Tests de validation pour src/thermodynamics/water_eos.jl
# Valeurs de référence :
#   - Kell (1975)                        → ρ(T, 1 bar)
#   - Helgeson (1969) / PHREEQC          → ε(T)
#   - Tanger & Helgeson (1988), Table 1  → Y, Q à 25 °C / 1 bar

@testset "water_density" begin
    ρ_25 = water_density(298.15, 1e5)
    @test isapprox(ρ_25, 0.9971; atol=1e-3)   # Kell (1975) : 0.99705 g/cm³

    ρ_100 = water_density(373.15, 1e5)
    @test isapprox(ρ_100, 0.9584; atol=5e-3)  # NIST : 0.9584 g/cm³

    # Monotonie : ρ décroît avec T à P = 1 bar
    @test ρ_100 < ρ_25

    # Compressibilité : ρ croît avec P
    ρ_25_500bar = water_density(298.15, 500e5)
    @test ρ_25_500bar > ρ_25

    # Ordre de grandeur physique
    @test 0.9 < ρ_25 < 1.1
end

@testset "water_dielectric" begin
    ε_25 = water_dielectric(298.15, 1e5)
    @test isapprox(ε_25, 78.54; atol=0.1)   # Helgeson (1969) donne exactement 78.54

    ε_100 = water_dielectric(373.15, 1e5)
    @test isapprox(ε_100, 55.7; atol=2.0)   # NIST : 55.7 à 100 °C, 1 bar

    # ε décroît avec T
    @test ε_100 < ε_25

    # ε croît légèrement avec P (électrostriction)
    ε_25_500bar = water_dielectric(298.15, 500e5)
    @test ε_25_500bar > ε_25

    @test ε_25 > 1.0
end

@testset "born_functions — signes et ordres de grandeur" begin
    bf = born_functions(298.15, 1e5)

    # ε et ρ cohérents
    @test isapprox(bf.ε, 78.54; atol=0.1)
    @test isapprox(bf.ρ, 0.9971; atol=1e-3)

    # Y = (1/ε²)(∂ε/∂T) < 0  (ε décroît avec T)
    # Tanger & Helgeson (1988), Table 1 : Y ≈ −5.81×10⁻⁵ K⁻¹
    @test bf.Y < 0.0
    @test isapprox(bf.Y, -5.81e-5; rtol=0.05)

    # Q = −(1/ε²)(∂ε/∂P) < 0  (ε croît légèrement avec P)
    # Tanger & Helgeson (1988), Table 1 : Q ≈ −3.57×10⁻¹⁵ Pa⁻¹
    @test bf.Q < 0.0
    @test isapprox(bf.Q, -3.57e-15; rtol=0.10)

    # X = ∂Y/∂T (fini, non nul)
    @test isfinite(bf.X)
    @test bf.X != 0.0
end

@testset "born_functions — variation avec T" begin
    bf_25  = born_functions(298.15, 1e5)
    bf_100 = born_functions(373.15, 1e5)

    @test bf_25.ε > bf_100.ε          # ε décroît avec T
    @test bf_25.Y != bf_100.Y          # Y varie avec T
    @test abs(bf_100.Y) > abs(bf_25.Y) # |Y| croît avec T (ε varie plus vite)
end
