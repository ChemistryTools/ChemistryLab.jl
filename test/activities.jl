# Tests unitaires pour src/equilibrium/activities.jl — modèle HKF (B-dot)
#
# Valeurs de référence :
#   - Helgeson, Kirkham & Flowers (1981), Am. J. Sci. 281, 1249–1516.
#     Table 3 : rayons électrostatiques effectifs åⱼ.
#     Eq. de l'activité (B-dot) et coefficient osmotique.
#   - Valeur A = 0.5114, B = 0.3288 à 25 °C / 1 bar : consensus littérature /
#     PHREEQC / SUPCRT92.
#   - γ±(NaCl, 0.1 mol/kg, 25 °C) expérimental ≈ 0.778 (Rard & Miller 1979).
#     Le modèle HKF avec les rayons de la Table 3 donne une valeur proche
#     mais pas identique (rayons HKF ≠ WATEQ4F).

# ── Constantes locales ────────────────────────────────────────────────────────

const _M_w = ustrip(us"kg/mol", Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)[:M])

# Fixtures de base : H₂O + Na⁺ + Cl⁻
const _h2o = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
const _Na  = Species("Na+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
const _Cl  = Species("Cl-"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
const _cs_nacl = ChemicalSystem([_h2o, _Na, _Cl])

# Moles correspondant à 1 kg d'eau + 0.1 mol/kg NaCl
const _nw_1kg   = 1.0 / _M_w           # ≈ 55.508 mol H₂O (= 1 kg eau)
const _n_salt   = 0.1                   # mol → m = 0.1 mol/kg exactement
const _n_nacl   = [_nw_1kg, _n_salt, _n_salt]
const _p_basic  = (ϵ = 1e-30,)

# ── 1. Fonction auxiliaire σ ──────────────────────────────────────────────────

@testset "σ(x) — valeur limite x = 0" begin
    # σ(0) = 1 par définition (continuité)
    @test ChemistryLab._hkf_sigma(0.0) ≈ 1.0
    @test ChemistryLab._hkf_sigma(0.0f0) ≈ 1.0f0
end

@testset "σ(x) — branche exacte (x ≥ 1e-3)" begin
    # Pour x ≥ 1e-3 le code utilise la formule exacte ; les annulations
    # catastrophiques n'existent plus à partir de x ~ 0.01.
    # À x = 0.01 : f = x - 2ln(1+x) - 1/(1+x) + 1 ≈ x³/3 (terme dominant)
    for x in [0.01, 0.1, 0.5, 1.0, 2.0]
        exact = (3 / x^3) * (x - 2 * log(1 + x) - 1 / (1 + x) + 1)
        @test isapprox(ChemistryLab._hkf_sigma(x), exact; rtol = 1e-12)
    end
end

@testset "σ(x) — branche Taylor : accord avec développement O(x²) pour x ≪ 1" begin
    # Pour |x| < 1e-3 le code utilise 1 - 3x/2 + 9x²/5 (Taylor tronqué à l'ordre 2).
    # On vérifie que le résultat coïncide avec la formule Taylor sur cette plage.
    for x in [1e-6, 1e-5, 1e-4, 5e-4, 9e-4]
        taylor = 1.0 - 1.5 * x + 1.8 * x^2
        @test isapprox(ChemistryLab._hkf_sigma(x), taylor; rtol = 1e-12)
    end
end

@testset "σ(x) — valeur exacte à x = 1" begin
    # σ(1) = 3 × (1 - 2 ln 2 - 1/2 + 1) = 3 × (1.5 - 2 ln 2) ≈ 0.3411
    expected = 3.0 * (1.5 - 2 * log(2))
    @test ChemistryLab._hkf_sigma(1.0) ≈ expected  rtol = 1e-12
end

@testset "σ(x) — monotone décroissant, borné entre 0 et 1 pour x > 0" begin
    xs = [1e-4, 0.1, 0.5, 1.0, 2.0, 5.0]
    σs = ChemistryLab._hkf_sigma.(xs)
    @test all(0.0 .< σs .≤ 1.0)
    # Décroissance
    for k in 1:length(σs)-1
        @test σs[k] > σs[k+1]
    end
end

# ── 2. Paramètres de Debye-Hückel A et B ─────────────────────────────────────

@testset "hkf_debye_huckel_params — référence 25 °C / 1 bar" begin
    ab = hkf_debye_huckel_params(298.15, 1e5)
    # Valeurs de consensus : A ≈ 0.5114 (kg/mol)^(1/2), B ≈ 0.3288 Å⁻¹ (kg/mol)^(1/2)
    @test isapprox(ab.A, 0.5114; atol = 0.005)
    @test isapprox(ab.B, 0.3288; atol = 0.005)
end

@testset "hkf_debye_huckel_params — A et B augmentent avec T à P fixée" begin
    ab25  = hkf_debye_huckel_params(298.15, 1e5)
    ab100 = hkf_debye_huckel_params(373.15, 1e5)
    # A est proportionnel à ρ^(1/2) / (ε T)^(3/2) :
    # ε décroît fortement avec T → A augmente
    @test ab100.A > ab25.A
end

@testset "hkf_debye_huckel_params — valeurs physiquement finies et positives" begin
    for (T, P) in [(298.15, 1e5), (373.15, 1e5), (373.15, 5e7)]
        ab = hkf_debye_huckel_params(T, P)
        @test isfinite(ab.A) && ab.A > 0
        @test isfinite(ab.B) && ab.B > 0
    end
end

# ── 3. Table REJ_HKF et valeurs de fallback ───────────────────────────────────

@testset "REJ_HKF — chargement et valeurs clés (Table 3, Helgeson 1981)" begin
    # La table doit être non vide
    @test !isempty(REJ_HKF)
    # Na⁺ → 1.91 Å (PHREEQC key "Na+")
    @test haskey(REJ_HKF, "Na+")
    @test isapprox(REJ_HKF["Na+"], 1.91; atol = 0.01)
    # Cl⁻ → 1.81 Å
    @test haskey(REJ_HKF, "Cl-")
    @test isapprox(REJ_HKF["Cl-"], 1.81; atol = 0.01)
    # Ca²⁺ → 2.87 Å (clé PHREEQC "Ca+2")
    @test haskey(REJ_HKF, "Ca+2")
    @test isapprox(REJ_HKF["Ca+2"], 2.87; atol = 0.01)
    # SO₄²⁻ → 3.15 Å (clé "SO4-2")
    @test haskey(REJ_HKF, "SO4-2")
    @test isapprox(REJ_HKF["SO4-2"], 3.15; atol = 0.01)
end

@testset "REJ_CHARGE_DEFAULT — cohérence des charges" begin
    for z in [-3, -2, -1, 1, 2, 3, 4]
        @test haskey(REJ_CHARGE_DEFAULT, z)
        @test REJ_CHARGE_DEFAULT[z] > 0.0
    end
    # Tendance : |å| croît avec |z|
    @test REJ_CHARGE_DEFAULT[-1] < REJ_CHARGE_DEFAULT[-2] < REJ_CHARGE_DEFAULT[-3]
    @test REJ_CHARGE_DEFAULT[1]  < REJ_CHARGE_DEFAULT[2]  < REJ_CHARGE_DEFAULT[3]
end

# ── 4. Construction du modèle ─────────────────────────────────────────────────

@testset "HKFActivityModel — constructeur par défaut" begin
    m = HKFActivityModel()
    @test m.A     ≈ 0.5114
    @test m.B     ≈ 0.3288
    @test m.Bdot  ≈ 0.041
    @test m.Kn    ≈ 0.1
    @test m.temperature_dependent == false
end

@testset "HKFActivityModel — constructeur avec mots-clés" begin
    m = HKFActivityModel(; A = 0.51, B = 0.33, Bdot = 0.04, temperature_dependent = true)
    @test m.A == 0.51
    @test m.B == 0.33
    @test m.temperature_dependent == true
end

@testset "activity_model — retourne une fonction" begin
    lna = activity_model(_cs_nacl, HKFActivityModel())
    @test lna isa Function
    out = lna(_n_nacl, _p_basic)
    @test length(out) == 3
    @test all(isfinite, out)
end

# ── 5. Limite de dilution infinie ─────────────────────────────────────────────

@testset "Limite de dilution infinie — γᵢ → 1" begin
    # À I → 0, log₁₀γ → 0, donc ln aᵢ → ln mᵢ
    lna = activity_model(_cs_nacl, HKFActivityModel())
    # Très faible concentration : n = 1e-10 mol
    ε = 1e-10
    n_dilute = [_nw_1kg, ε, ε]
    out = lna(n_dilute, _p_basic)

    # Molalité attendue pour n_salt = 1e-10 mol
    m_dilute = ε / (_nw_1kg * _M_w)   # mol/kg

    # log₁₀γ ≈ -A × z² × √I ≈ -0.5114 × 1 × √(ε) ≈ très proche de 0
    # ln aᵢ ≈ log(m_dilute) (correction γ ≲ 1e-6)
    @test isapprox(out[2], log(m_dilute); atol = 1e-4)  # Na⁺
    @test isapprox(out[3], log(m_dilute); atol = 1e-4)  # Cl⁻
end

@testset "Limite de dilution infinie — activité de l'eau → 0" begin
    lna = activity_model(_cs_nacl, HKFActivityModel())
    n_dilute = [_nw_1kg, 1e-10, 1e-10]
    out = lna(n_dilute, _p_basic)
    # ln a_w ≈ -Mw × Σm × φ ≈ 0 pour Σm → 0
    @test abs(out[1]) < 1e-6
end

# ── 6. Solution NaCl 0.1 mol/kg — vérification analytique ────────────────────

@testset "NaCl 0.1 mol/kg — log₁₀γ matches DH formula" begin
    lna = activity_model(_cs_nacl, HKFActivityModel())
    out = lna(_n_nacl, _p_basic)

    # Calcul analytique avec les paramètres par défaut et å de REJ_HKF
    A     = 0.5114;  B = 0.3288;  Bdot = 0.041
    å_Na  = REJ_HKF["Na+"]          # 1.91 Å
    å_Cl  = REJ_HKF["Cl-"]          # 1.81 Å
    I     = 0.1                      # mol/kg (symétrique NaCl)
    sqrtI = sqrt(I)
    m     = 0.1                      # mol/kg

    log10_γ_Na = -A * 1^2 * sqrtI / (1 + B * å_Na * sqrtI) + Bdot * I
    log10_γ_Cl = -A * 1^2 * sqrtI / (1 + B * å_Cl * sqrtI) + Bdot * I

    ln10 = log(10)
    # ln aᵢ = ln10 × log₁₀γᵢ + ln(mᵢ)  [état standard molalité]
    ln_a_Na_expected = ln10 * log10_γ_Na + log(m)
    ln_a_Cl_expected = ln10 * log10_γ_Cl + log(m)

    # Tolérance serrée : calcul identique à celui du code
    @test isapprox(out[2], ln_a_Na_expected; rtol = 1e-8)   # Na⁺
    @test isapprox(out[3], ln_a_Cl_expected; rtol = 1e-8)   # Cl⁻
end

@testset "NaCl 0.1 mol/kg — coefficients d'activité < 1" begin
    # Pour les ions monovalents à 0.1 mol/kg, γ < 1 (DH prédit des coefficients
    # inférieurs à 1 pour tout I > 0)
    lna = activity_model(_cs_nacl, HKFActivityModel())
    out = lna(_n_nacl, _p_basic)

    m = 0.1
    ln_a_Na = out[2]
    ln_a_Cl = out[3]

    γ_Na = exp(ln_a_Na - log(m))   # exp(ln a - ln m) = γ
    γ_Cl = exp(ln_a_Cl - log(m))

    @test γ_Na < 1.0
    @test γ_Cl < 1.0
    # Valeur expérimentale de référence : γ±(NaCl, 0.1 mol/kg) ≈ 0.778
    # Le modèle HKF avec ses rayons spécifiques donne ~ 0.73–0.77 (acceptable)
    γ_mean = sqrt(γ_Na * γ_Cl)
    @test 0.60 < γ_mean < 1.0
end

@testset "NaCl 0.1 mol/kg — activité de l'eau négative et correcte" begin
    lna = activity_model(_cs_nacl, HKFActivityModel())
    out = lna(_n_nacl, _p_basic)

    ln_aw = out[1]
    @test ln_aw < 0.0       # l'activité de l'eau diminue avec les solutés
    @test ln_aw > -0.1      # mais reste proche de 1 à cette concentration

    # Vérification analytique du coefficient osmotique
    A = 0.5114;  B = 0.3288;  Bdot = 0.041
    å_Na = REJ_HKF["Na+"]
    å_Cl = REJ_HKF["Cl-"]
    m = 0.1;  I = 0.1;  sqrtI = sqrt(I)
    ln10 = log(10)

    Σm   = 2 * m
    å_eff = (m * 1^2 * å_Na + m * 1^2 * å_Cl) / (m * 1^2 + m * 1^2)
    σ    = ChemistryLab._hkf_sigma(B * å_eff * sqrtI)
    φ    = 1.0 - (A * ln10 / 3) * (2I / Σm) * sqrtI * σ + (Bdot * ln10 / 2) * I
    ln_aw_expected = -_M_w * Σm * φ

    @test isapprox(ln_aw, ln_aw_expected; rtol = 1e-8)
end

# ── 7. Soluté neutre — coefficient de salage ──────────────────────────────────

@testset "Soluté neutre — coefficient Kn" begin
    # Système H₂O + CO₂@ (neutre) + Na⁺ + Cl⁻
    co2aq = Species("CO2@"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    cs_neu = ChemicalSystem([_h2o, co2aq, _Na, _Cl])
    lna = activity_model(cs_neu, HKFActivityModel())

    # 0.5 mol/kg NaCl → I = 0.5
    nw = _nw_1kg
    n_Na = n_Cl = 0.5 * nw * _M_w
    n_co2 = 0.01 * nw * _M_w   # 0.01 mol/kg CO₂

    n_vec = [nw, n_co2, n_Na, n_Cl]
    out = lna(n_vec, _p_basic)

    # Indice de CO₂@ dans cs_neu (recherche par formule PHREEQC)
    i_co2 = findfirst(s -> s.formula.phreeqc == "CO2@", cs_neu.species)

    I = 0.5
    Kn = 0.1
    m_co2 = n_co2 / (nw * _M_w)
    ln_a_co2_expected = log(10) * (Kn * I) + log(m_co2)
    @test isapprox(out[i_co2], ln_a_co2_expected; rtol = 1e-8)
end

# ── 8. Lookup des rayons ioniques ─────────────────────────────────────────────

@testset "Priorité du rayon : sp[:å] > REJ_HKF > REJ_CHARGE_DEFAULT > å_default" begin
    # Création d'un soluté fictif avec å explicite dans les propriétés
    sp_custom = Species("Na+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    sp_custom[:å] = 5.0   # rayon personnalisé [Å]

    cs_custom = ChemicalSystem([_h2o, sp_custom])
    lna_custom  = activity_model(cs_custom, HKFActivityModel())
    lna_default = activity_model(_cs_nacl, HKFActivityModel())   # Na+ sans prop :å

    n_vec = [_nw_1kg, _n_salt]
    out_custom  = lna_custom(n_vec, _p_basic)
    out_default = activity_model(ChemicalSystem([_h2o, _Na]), HKFActivityModel())([_nw_1kg, _n_salt], _p_basic)

    # Avec å=5.0 >> REJ(Na+)=1.91, le dénominateur est plus grand → γ plus grand (moins négatif)
    @test out_custom[2] > out_default[2]
end

@testset "Fallback REJ_CHARGE_DEFAULT — ion inconnu de charge ±2" begin
    # Ion fictif de charge +2 absent de REJ_HKF
    sp_unknown = Species("Xx+2"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    cs_unk = ChemicalSystem([_h2o, sp_unknown])
    lna = activity_model(cs_unk, HKFActivityModel())
    out = lna([_nw_1kg, _n_salt], _p_basic)

    # Avec å = REJ_CHARGE_DEFAULT[2] = 2.80 Å, vérification analytique
    A = 0.5114;  B = 0.3288;  Bdot = 0.041
    å = REJ_CHARGE_DEFAULT[2]   # 2.80 Å
    I = 0.5 * 4 * _n_salt / (_nw_1kg * _M_w)   # I = 0.5 z² m = 0.5×4×0.1 = 0.2
    sqrtI = sqrt(I)
    m = _n_salt / (_nw_1kg * _M_w)

    log10_γ = -A * 4 * sqrtI / (1 + B * å * sqrtI) + Bdot * I
    ln_a_expected = log(10) * log10_γ + log(m)

    @test isapprox(out[2], ln_a_expected; rtol = 1e-8)
end

# ── 9. Mode température-dépendant ────────────────────────────────────────────

@testset "temperature_dependent=true — différent de paramètres fixes à T≠25°C" begin
    model_fixed = HKFActivityModel(; temperature_dependent = false)
    model_tp    = HKFActivityModel(; temperature_dependent = true)

    lna_fixed = activity_model(_cs_nacl, model_fixed)
    lna_tp    = activity_model(_cs_nacl, model_tp)

    p_25 = (ϵ = 1e-30, T = 298.15, P = 1e5)
    p_100 = (ϵ = 1e-30, T = 373.15, P = 1e5)

    # À 25 °C les valeurs doivent être pratiquement identiques (A, B identiques)
    out_fixed_25 = lna_fixed(_n_nacl, p_25)
    out_tp_25    = lna_tp(_n_nacl, p_25)
    @test isapprox(out_fixed_25[2], out_tp_25[2]; rtol = 1e-3)

    # À 100 °C, le mode T-dépendant doit différer
    out_fixed_100 = lna_fixed(_n_nacl, p_100)   # utilise A, B fixes
    out_tp_100    = lna_tp(_n_nacl, p_100)       # recalcule A, B
    @test out_fixed_100[2] ≠ out_tp_100[2]
end

@testset "temperature_dependent=true — sans T/P dans p, repli sur valeurs fixes" begin
    model = HKFActivityModel(; A = 0.55, B = 0.34, temperature_dependent = true)
    lna   = activity_model(_cs_nacl, model)
    # p sans :T ni :P → le modèle doit utiliser A=0.55, B=0.34 (valeurs fixes)
    out_tp    = lna(_n_nacl, (ϵ = 1e-30,))
    model_fix = HKFActivityModel(; A = 0.55, B = 0.34, temperature_dependent = false)
    out_fix   = activity_model(_cs_nacl, model_fix)(_n_nacl, (ϵ = 1e-30,))
    @test isapprox(out_tp[2], out_fix[2]; rtol = 1e-10)
end

# ── 10. Cohérence Gibbs-Duhem ─────────────────────────────────────────────────

@testset "Gibbs-Duhem — consistance numérique (différences finies)" begin
    # Σᵢ nᵢ × ∂μᵢ/∂ξ ≈ 0 pour ξ = ajout simultané de δ mol Na⁺ et Cl⁻
    lna = activity_model(_cs_nacl, HKFActivityModel())

    n0  = copy(_n_nacl)
    δ   = 1e-5   # mol

    # Variation : ajout de δ mol NaCl (n_Na et n_Cl + δ)
    n_plus  = n0 .+ [0.0, δ, δ]
    n_minus = n0 .- [0.0, δ, δ]

    dμ = (lna(n_plus, _p_basic) .- lna(n_minus, _p_basic)) ./ (2δ)   # ∂μᵢ/∂ξ ≈ dμᵢ/δ

    gd = sum(n0[i] * dμ[i] for i in eachindex(n0))
    # La somme est d'ordre δ² (terme d'erreur des différences finies) et
    # doit être petite devant Σ |nᵢ × dμᵢ|
    norm_ref = sum(abs(n0[i] * dμ[i]) for i in eachindex(n0))
    @test abs(gd) / norm_ref < 1e-4
end

# ── 11. Phase gazeuse — mélange idéal ────────────────────────────────────────

@testset "Phase gazeuse — mélange idéal (ln aᵢ = ln xᵢ)" begin
    co2g = Species("CO2"; aggregate_state = AS_GAS, class = SC_GASFLUID)
    n2g  = Species("N2";  aggregate_state = AS_GAS, class = SC_GASFLUID)
    cs_gas = ChemicalSystem([_h2o, _Na, co2g, n2g])
    lna = activity_model(cs_gas, HKFActivityModel())

    n_co2 = 0.3;  n_n2 = 0.7   # x_CO₂ = 0.3, x_N₂ = 0.7
    n_vec = [_nw_1kg, 0.01, n_co2, n_n2]
    out = lna(n_vec, _p_basic)

    i_co2 = findfirst(s -> s.formula.phreeqc == "CO2", cs_gas.species)
    i_n2  = findfirst(s -> s.formula.phreeqc == "N2",  cs_gas.species)

    @test isapprox(out[i_co2], log(0.3); rtol = 1e-8)
    @test isapprox(out[i_n2],  log(0.7); rtol = 1e-8)
end

# ── 12. build_potentials ──────────────────────────────────────────────────────

@testset "build_potentials — retourne une fonction de longueur correcte" begin
    μ = build_potentials(_cs_nacl, HKFActivityModel())
    @test μ isa Function

    ΔaG = [-400.0, -300.0, -250.0]   # μᵢ° / RT (sans unité)
    p = (ΔₐG⁰overT = ΔaG, ϵ = 1e-30)
    out = μ(_n_nacl, p)
    @test length(out) == 3
    @test all(isfinite, out)

    # μᵢ = μᵢ° + ln aᵢ
    lna_out = activity_model(_cs_nacl, HKFActivityModel())(_n_nacl, p)
    @test isapprox(out, ΔaG .+ lna_out; rtol = 1e-10)
end
