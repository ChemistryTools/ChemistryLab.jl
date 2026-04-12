using ForwardDiff

# ── arrhenius_rate_constant ───────────────────────────────────────────────────

@testset "arrhenius_rate_constant" begin

    k = arrhenius_rate_constant(5.0e-4, 62000.0)

    @test k isa NumericFunc

    # At T_ref = 298.15, k(T_ref) = k₀
    @test isapprox(k(; T = 298.15), 5.0e-4; rtol = 1.0e-10)

    # Arrhenius: higher T → higher k
    @test k(; T = 350.0) > k(; T = 298.15)
    @test k(; T = 250.0) < k(; T = 298.15)

    # AD smoke-test: dk/dT should be finite and positive
    dkdT = ForwardDiff.derivative(T -> arrhenius_rate_constant(5.0e-4, 62000.0)(; T = T), 298.15)
    @test isfinite(dkdT)
    @test dkdT > 0

    # AD w.r.t. Ea: at T > T_ref, larger Ea increases k
    # dk/dEa = k(T) × (−1/R) × (1/T − 1/T_ref); T=350 > T_ref=298.15 → (1/T−1/T_ref)<0 → dk/dEa>0
    dk_dEa = ForwardDiff.derivative(Ea -> arrhenius_rate_constant(5.0e-4, Ea)(; T = 350.0), 62000.0)
    @test isfinite(dk_dEa)
    @test dk_dEa > 0   # at T > T_ref the exponent is positive; larger Ea → larger k

    # AD w.r.t. k₀
    dk_dk0 = ForwardDiff.derivative(k₀ -> arrhenius_rate_constant(k₀, 62000.0)(; T = 298.15), 5.0e-4)
    @test isapprox(dk_dk0, 1.0; rtol = 1.0e-10)  # k(T_ref) = k₀ linearly

end

# ── KINETICS_RATE_FACTORIES ───────────────────────────────────────────────────

@testset "KINETICS_RATE_FACTORIES" begin

    @test haskey(KINETICS_RATE_FACTORIES, :arrhenius)

    factory = KINETICS_RATE_FACTORIES[:arrhenius]
    k = factory(; k₀ = 1.0e-6, Ea = 40000.0, T_ref = 298.15, R_gas = 8.31446)
    @test k isa AbstractFunc
    @test isapprox(k(; T = 298.15), 1.0e-6; rtol = 1.0e-10)

end

# ── saturation_ratio ──────────────────────────────────────────────────────────

@testset "saturation_ratio" begin

    # At equilibrium (IAP = K): Ω = 1
    # For reaction A ⇌ B:  stoich = [-1, 1],
    #   ln IAP = ln a_B - ln a_A,  ln K = -(ΔₐG⁰_B - ΔₐG⁰_A)/RT
    # Set lna_A = lna_B = 0 and ΔₐG⁰/RT both zero → Ω = exp(0 - 0) = 1
    stoich = [-1.0, 1.0]
    lna = [0.0, 0.0]
    G_over_T = [0.0, 0.0]
    Ω = saturation_ratio(stoich, lna, G_over_T)
    @test isapprox(Ω, 1.0; rtol = 1.0e-12)

    # Ω < 1: system undersaturated w.r.t. reaction forward direction
    # lna_B < lna_A (less product than equilibrium requires)
    lna2 = [-1.0, 0.0]    # a_A = exp(-1), a_B = 1
    Ω2 = saturation_ratio(stoich, lna2, G_over_T)
    @test Ω2 > 1.0   # excess A → forward reaction is favoured

    # AD smoke-test
    dΩ = ForwardDiff.derivative(
        x -> saturation_ratio(stoich, [x, 0.0], G_over_T),
        0.0,
    )
    @test isfinite(dΩ)

end

# ── TransitionStateRateModel ──────────────────────────────────────────────────

@testset "TransitionStateRateModel" begin

    k_neutral = arrhenius_rate_constant(1.549e-6, 23500.0)
    k_acid = arrhenius_rate_constant(5.012e-1, 14400.0)

    mech_neutral = RateMechanism(k_neutral, 1.0, 1.0)
    mech_acid = RateMechanism(k_acid, 1.0, 1.0, [RateModelCatalyst("H+", 1.0)])

    model = TransitionStateRateModel([mech_neutral, mech_acid])
    @test model isa TransitionStateRateModel

    # At equilibrium (Ω = 1): rate ≈ 0
    r_eq = model(; T = 298.15, Ω = 1.0, A_surface = 1.0)
    @test abs(r_eq) < 1.0e-12

    # Dissolution (Ω < 1): positive rate
    r_diss = model(; T = 298.15, Ω = 0.0, A_surface = 1.0)
    @test r_diss > 0

    # Precipitation (Ω > 1): negative rate
    r_prec = model(; T = 298.15, Ω = 2.0, A_surface = 1.0)
    @test r_prec < 0

    # Rate scales with surface area
    r1 = model(; T = 298.15, Ω = 0.5, A_surface = 1.0)
    r2 = model(; T = 298.15, Ω = 0.5, A_surface = 2.0)
    @test isapprox(r2, 2 * r1; rtol = 1.0e-10)

    # AD smoke-test w.r.t. Ω
    drdΩ = ForwardDiff.derivative(Ω -> model(; T = 298.15, Ω = Ω, A_surface = 1.0), 0.5)
    @test isfinite(drdΩ)
    @test drdΩ < 0   # higher saturation → lower dissolution rate

    # AD w.r.t. T (via closure over k)
    drdT = ForwardDiff.derivative(T -> model(; T = T, Ω = 0.5, A_surface = 1.0), 298.15)
    @test isfinite(drdT)
    @test drdT > 0   # higher T → faster dissolution

end

# ── FirstOrderRateModel ───────────────────────────────────────────────────────

@testset "FirstOrderRateModel" begin

    k = arrhenius_rate_constant(1.0e-7, 40000.0)
    rm = FirstOrderRateModel(k)

    @test rm(; T = 298.15, Ω = 1.0, A_surface = 1.0) ≈ 0  atol = 1.0e-12
    @test rm(; T = 298.15, Ω = 0.5, A_surface = 1.0) > 0
    # Extra kwargs ignored
    @test rm(; T = 298.15, Ω = 0.5, A_surface = 1.0, lna_dict = nothing) ==
        rm(; T = 298.15, Ω = 0.5, A_surface = 1.0)

    drdΩ = ForwardDiff.derivative(Ω -> rm(; T = 298.15, Ω = Ω, A_surface = 1.0), 0.5)
    @test isfinite(drdΩ)

end

# ── Surface area models ───────────────────────────────────────────────────────

@testset "SurfaceArea models" begin

    fixed = FixedSurfaceArea(0.5)
    @test surface_area(fixed, 0.01, 0.1) ≈ 0.5
    @test surface_area(fixed, 0.0, 0.1) ≈ 0.5   # constant regardless of n

    bet = BETSurfaceArea(90.0)   # 90 m²/kg
    @test surface_area(bet, 0.01, 0.1) ≈ 90.0 * 0.01 * 0.1   # = 0.09 m²
    @test surface_area(bet, 0.0, 0.1) ≈ 0.0   # zero moles → zero area

    # AD through BET surface area
    dA = ForwardDiff.derivative(n -> surface_area(bet, n, 0.1), 0.01)
    @test isfinite(dA)
    @test dA > 0

end

# ── ParrotKillohRateModel ────────────────────────────────────────────────────

@testset "ParrotKillohRateModel" begin

    # ── Construction ─────────────────────────────────────────────────────────
    rm = PK_PARAMS_C3S
    @test rm isa ParrotKillohRateModel
    @test rm.K₁ ≈ 1.5
    @test rm.T_ref ≈ 293.15
    @test rm.α_max ≈ 1.0

    # Custom α_max
    rm2 = ParrotKillohRateModel(
        1.5, 3.3, 0.018, 2.5, 0.0024, 4.0, 0.5, 41_570.0;
        T_ref = 293.15, α_max = 0.85,
    )
    @test rm2.α_max ≈ 0.85

    # ── Positivity and monotonicity ───────────────────────────────────────────
    # Rate should be positive (dissolution) for any α ∈ (0, α_max)
    T_K = 293.15   # at reference temperature → Aₜ = 1
    n0 = 1.0      # mol
    r_early = rm(; T = T_K, n_current = 0.99 * n0, n_initial = n0)   # α ≈ 0.01
    r_middle = rm(; T = T_K, n_current = 0.5 * n0, n_initial = n0)   # α = 0.50
    r_late = rm(; T = T_K, n_current = 0.01 * n0, n_initial = n0)   # α ≈ 0.99
    @test r_early > 0
    @test r_middle > 0
    @test r_late > 0

    # At α = α_max (mineral fully consumed), rate should be negligible
    r_max = rm(; T = T_K, n_current = 0.0, n_initial = n0)
    @test r_max < r_early   # rate should be small (clamped at α_max - ε)

    # ── Temperature dependence ────────────────────────────────────────────────
    # Higher T → higher rate (Arrhenius)
    r_low_T = rm(; T = 273.15, n_current = 0.5 * n0, n_initial = n0)
    r_high_T = rm(; T = 323.15, n_current = 0.5 * n0, n_initial = n0)
    @test r_high_T > r_low_T

    # At T_ref: Aₜ = 1 → check consistency with manual computation
    # For C3S at α=0 (ξ=0): r_NG = K₁/N₁ / 86400, r_I = K₂/86400, r_D → ∞
    # → min(max(r_NG, r_I), r_D) = max(r_NG, r_I) = r_NG (since K₁/N₁ > K₂)
    # rate at n_current ≈ n0 is ≈ n0 * K₁/N₁ / 86400 / (1 + 0) = n0 * 1.5/3.3 / 86400
    expected_early = n0 * (rm.K₁ / rm.N₁) / 86400  # ξ≈0 → NG dominates
    @test isapprox(r_early, expected_early; rtol = 0.05)  # 5% tol (ξ=0.01 is not exactly 0)

    # ── Predefined constants for all 4 clinker phases ────────────────────────
    @test PK_PARAMS_C2S isa ParrotKillohRateModel
    @test PK_PARAMS_C3A isa ParrotKillohRateModel
    @test PK_PARAMS_C4AF isa ParrotKillohRateModel
    for rm_k in (PK_PARAMS_C3S, PK_PARAMS_C2S, PK_PARAMS_C3A, PK_PARAMS_C4AF)
        r = rm_k(; T = 293.15, n_current = 0.5, n_initial = 1.0)
        @test r > 0
        @test isfinite(r)
    end

    # ── Interface compatibility: extra kwargs silently ignored ────────────────
    r_with_extra = rm(;
        T = T_K, n_current = 0.5 * n0, n_initial = n0,
        Ω = 0.5, A_surface = 1.0, lna_dict = nothing,
    )
    r_bare = rm(; T = T_K, n_current = 0.5 * n0, n_initial = n0)
    @test r_with_extra ≈ r_bare

    # ── AD smoke-tests ────────────────────────────────────────────────────────
    # w.r.t. n_current
    dr_dnc = ForwardDiff.derivative(
        nc -> PK_PARAMS_C3S(; T = 293.15, n_current = nc, n_initial = 1.0),
        0.5,
    )
    @test isfinite(dr_dnc)

    # w.r.t. T
    dr_dT = ForwardDiff.derivative(
        T -> PK_PARAMS_C3S(; T = T, n_current = 0.5, n_initial = 1.0),
        293.15,
    )
    @test isfinite(dr_dT)
    @test dr_dT > 0   # higher T → faster rate

    # w.r.t. K₁ (construct fresh model to propagate Dual)
    dr_dK₁ = ForwardDiff.derivative(0.0) do K₁_val
        m = ParrotKillohRateModel(
            K₁_val, 3.3, 0.018, 2.5, 0.0024, 4.0, 0.5, 41_570.0;
            T_ref = 293.15, α_max = 1.0,
        )
        m(; T = 293.15, n_current = 0.5, n_initial = 1.0)
    end  # K₁=0 → rate=0; derivative>0
    @test isfinite(dr_dK₁)
    @test dr_dK₁ ≥ 0

end

# ── KineticReaction convenience constructor ────────────────────────────────────

@testset "KineticReaction convenience constructor" begin

    # ── Build a minimal ChemicalSystem with one mineral + solvent ─────────────
    calcite = Species("Calcite"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    h2o = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    ca2p = Species("Ca+2"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)

    cs = ChemicalSystem([calcite, h2o, ca2p])
    n_sp = length(cs.species)   # = 3

    k_neutral = arrhenius_rate_constant(1.549e-6, 23500.0)
    rm = FirstOrderRateModel(k_neutral)
    sm = FixedSurfaceArea(0.5)

    # ── Constructor finds species by phreeqc formula ───────────────────────────
    kr = KineticReaction(cs, "Calcite", rm, sm)
    @test kr isa KineticReaction
    @test kr.idx_mineral == 1                            # Calcite is species[1]
    @test length(kr.stoich) == n_sp
    @test kr.stoich[1] ≈ -1.0                            # mineral
    @test all(iszero, kr.stoich[2:end])                  # other species = 0
    @test kr.rate_model === rm
    @test kr.surface_model === sm

    # ── molar_mass dispatch on AbstractSpecies path ───────────────────────────
    # Calcite has no :M set → falls back to 0.1; after setting it, reads correctly
    calcite2 = Species("CaCO3"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    calcite2.properties[:M] = 0.10009u"kg/mol"
    cs2 = ChemicalSystem([calcite2, h2o])
    kr2 = KineticReaction(cs2, "CaCO3", rm, sm)
    @test isapprox(molar_mass(kr2), 0.10009; rtol = 1.0e-6)

    # ── Custom stoichiometry keyword ───────────────────────────────────────────
    stoich_custom = [-1.0, 0.0, 1.0]
    kr_custom = KineticReaction(cs, "Calcite", rm, sm; stoich = stoich_custom)
    @test kr_custom.stoich ≈ stoich_custom

    # ── DimensionMismatch for wrong-length stoich ──────────────────────────────
    @test_throws DimensionMismatch KineticReaction(cs, "Calcite", rm, sm; stoich = [-1.0])

    # ── ArgumentError for unknown species ─────────────────────────────────────
    @test_throws ArgumentError KineticReaction(cs, "Quartz", rm, sm)

    # ── Low-level constructor: explicit reaction, index, stoich ───────────────
    rxn = Reaction([calcite, ca2p]; symbol = "calcite dissolution")
    kr_low = KineticReaction(rxn, rm, sm, 1, [-1.0, 0.0, 1.0])
    @test kr_low.idx_mineral == 1
    @test kr_low.stoich[3] ≈ 1.0

    # Low-level: validate idx_mineral > 0
    @test_throws ArgumentError KineticReaction(rxn, rm, sm, 0, [-1.0, 0.0, 1.0])
    # Low-level: validate non-empty stoich
    @test_throws ArgumentError KineticReaction(rxn, rm, sm, 1, Float64[])

end
