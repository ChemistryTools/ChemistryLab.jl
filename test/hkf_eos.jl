# Tests de validation pour src/thermodynamics/hkf_eos.jl
# Valeurs de référence : SUPCRT92, Tanger & Helgeson (1988), Helgeson et al. (1981)
#
# Espèce de référence : Na⁺
# Paramètres SUPCRT92 (Helgeson et al. 1981, Table 2) :
#   a₁ = 0.1839 cal/mol/bar, a₂ = -228.59 cal/mol,
#   a₃ = 3.2564 cal·K/mol/bar, a₄ = -27884 cal·K/mol,
#   c₁ = 18.18 cal/mol/K, c₂ = -29810 cal·K/mol, ω = 33060 cal/mol
#   Gr = -62591 cal/mol, Hr = -57433 cal/mol, Sr = 13.96 cal/mol/K
#
# Valeurs SUPCRT92 (Table 5 de Tanger & Helgeson 1988 / Table B1 de Shock & Helgeson 1988)
# à 25 °C / 1 bar : G° = -261.9 kJ/mol, Cp° = 38.2 J/mol/K

# ── Construction des paramètres Na⁺ ──────────────────────────────────────────

# Paramètres SUPCRT92 en unités calorie/bar (Helgeson et al. 1981)
const _Na_params = HKFParams(;
    a₁ = 0.1839,
    a₂ = -228.59,
    a₃ = 3.2564,
    a₄ = -27884.0,
    c₁ = 18.18,
    c₂ = -29810.0,
    ω  = 33060.0,
    units = :supcrt,
)

const _Na_ref = HKFReferenceState(
    Gr = -62591.0 * 4.184,   # cal/mol → J/mol
    Hr = -57433.0 * 4.184,
    Sr = 13.96 * 4.184,      # cal/mol/K → J/mol/K
)

@testset "HKFParams — conversion d'unités" begin
    # a₁ : 0.1839 cal/mol/bar = 0.1839 × 4.184 / 1e5 J/mol/Pa
    @test isapprox(_Na_params.a₁, 0.1839 * 4.184 / 1e5; rtol = 1e-6)
    # ω : 33060 cal/mol → J/mol
    @test isapprox(_Na_params.ω, 33060.0 * 4.184; rtol = 1e-6)
end

@testset "HKFParams — constructeur :si" begin
    p_si = HKFParams(;
        a₁ = _Na_params.a₁,
        a₂ = _Na_params.a₂,
        a₃ = _Na_params.a₃,
        a₄ = _Na_params.a₄,
        c₁ = _Na_params.c₁,
        c₂ = _Na_params.c₂,
        ω  = _Na_params.ω,
        units = :si,
    )
    @test p_si.a₁ === _Na_params.a₁
    @test p_si.ω  === _Na_params.ω
end

@testset "hkf_G — point de référence Na⁺" begin
    # G° au point de référence doit être exactement Gr (erreurs numériques Born < 1 J/mol)
    G_ref = hkf_G(_Na_params, _Na_ref, HKF_Tr, HKF_Pr)
    @test isapprox(G_ref, _Na_ref.Gr; atol = 10.0)   # tolérance 10 J/mol

    # Valeur absolue : −261.9 kJ/mol (Shock & Helgeson 1988)
    @test isapprox(G_ref, -261.9e3; atol = 200.0)    # ±200 J/mol
end

@testset "hkf_H — point de référence Na⁺" begin
    H_ref = hkf_H(_Na_params, _Na_ref, HKF_Tr, HKF_Pr)
    @test isapprox(H_ref, _Na_ref.Hr; atol = 10.0)
end

@testset "hkf_S — point de référence Na⁺" begin
    S_ref = hkf_S(_Na_params, _Na_ref, HKF_Tr, HKF_Pr)
    @test isapprox(S_ref, _Na_ref.Sr; atol = 0.1)
end

@testset "hkf_Cp — point de référence Na⁺" begin
    # Cp° à 25°C/1 bar ≈ 38.2 J/mol/K (SUPCRT92 / Shock & Helgeson 1988).
    # La tolérance élargie (±5 J/mol/K) tient compte du fait que le modèle diélectrique
    # utilisé ici (Helgeson 1969) diffère légèrement du modèle HK1974 de SUPCRT92,
    # ce qui affecte le terme Born X dans Cp°.
    Cp_ref = hkf_Cp(_Na_params, HKF_Tr, HKF_Pr)
    @test isapprox(Cp_ref, 38.2; atol = 5.0)         # ±5 J/mol/K (modèle H69 vs HK1974)
end

@testset "hkf_V — signe et ordre de grandeur Na⁺" begin
    V_ref = hkf_V(_Na_params, HKF_Tr, HKF_Pr)
    # Volume molaire partiel de Na⁺ ≈ -1.2 cm³/mol = -1.2e-6 m³/mol (Shock & Helgeson 1988)
    # (volume partiel négatif possible pour ions très petits)
    @test isfinite(V_ref)
    @test abs(V_ref) < 1e-3    # < 1 L/mol (ordre de grandeur physique)
end

@testset "Consistance thermodynamique — relations de Maxwell" begin
    # ∂G/∂T ≈ −S  et  ∂G/∂P ≈ V  (différences finies centrées)
    ΔT = 1.0    # K
    ΔP = 1.0e4  # Pa

    for (T, P) in [(298.15, 1e5), (373.15, 1e5), (373.15, 5e7)]
        G₊T = hkf_G(_Na_params, _Na_ref, T + ΔT, P)
        G₋T = hkf_G(_Na_params, _Na_ref, T - ΔT, P)
        dGdT = (G₊T - G₋T) / (2 * ΔT)
        S = hkf_S(_Na_params, _Na_ref, T, P)
        @test isapprox(dGdT, -S; rtol = 1e-4)      # ∂G/∂T = −S

        G₊P = hkf_G(_Na_params, _Na_ref, T, P + ΔP)
        G₋P = hkf_G(_Na_params, _Na_ref, T, P - ΔP)
        dGdP = (G₊P - G₋P) / (2 * ΔP)
        V = hkf_V(_Na_params, T, P)
        @test isapprox(dGdP, V; rtol = 1e-4)       # ∂G/∂P = V
    end
end

@testset "Variation de G avec T" begin
    G_25  = hkf_G(_Na_params, _Na_ref, 298.15, 1e5)
    G_100 = hkf_G(_Na_params, _Na_ref, 373.15, 1e5)
    # G doit varier de manière finie et physiquement raisonnable
    @test isfinite(G_25)
    @test isfinite(G_100)
    @test abs(G_100 - G_25) < 50e3    # variation < 50 kJ/mol entre 25 et 100 °C
end

@testset "build_hkf_functions — interface" begin
    fns = build_hkf_functions(_Na_params, _Na_ref)

    @test haskey(fns, :ΔₐG⁰)
    @test haskey(fns, :ΔₐH⁰)
    @test haskey(fns, :S⁰)
    @test haskey(fns, :Cp⁰)
    @test haskey(fns, :V⁰)

    # Appel par mot-clé (interface compatible ThermoFunction)
    G = fns[:ΔₐG⁰](T = 298.15, P = 1e5)
    @test isapprox(G, -261.9e3; atol = 200.0)

    # Appel avec unit=true retourne une Quantity
    G_q = fns[:ΔₐG⁰](T = 298.15, P = 1e5, unit = true)
    @test G_q isa Quantity
end

@testset "HKFThermoFunction — propriété inconnue" begin
    f_bad = HKFThermoFunction(_Na_params, _Na_ref, :foo, DEFAULT_WATER_EOS)
    @test_throws ArgumentError f_bad(T = 298.15, P = 1e5)
end
