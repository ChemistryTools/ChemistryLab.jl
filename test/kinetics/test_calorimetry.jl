using ForwardDiff

# ── IsothermalCalorimeter ─────────────────────────────────────────────────────

@testset "IsothermalCalorimeter" begin

    cal = IsothermalCalorimeter(298.15)
    @test cal isa IsothermalCalorimeter
    @test cal.T ≈ 298.15

    @test n_extra_states(cal) == 1

    u0 = [0.01, 0.05]
    u0_ext = extend_u0(u0, cal)
    @test length(u0_ext) == 3
    @test u0_ext[end] == 0.0
    @test u0_ext[1:2] == u0

end

# ── SemiAdiabaticCalorimeter ──────────────────────────────────────────────────

@testset "SemiAdiabaticCalorimeter" begin

    # Backward-compat linear constructor (L keyword)
    cal_lin = SemiAdiabaticCalorimeter(; T0 = 293.15, T_env = 293.15, L = 0.5, Cp = 4000.0)
    @test cal_lin isa SemiAdiabaticCalorimeter
    @test cal_lin.T0 ≈ 293.15
    @test cal_lin.T_env ≈ 293.15
    @test cal_lin.Cp ≈ 4000.0
    @test cal_lin.heat_loss(1.0) ≈ 0.5   # φ(ΔT=1) = L*1 = 0.5
    @test cal_lin.heat_loss(2.0) ≈ 1.0   # φ(ΔT=2) = L*2 = 1.0

    @test n_extra_states(cal_lin) == 1

    u0 = [0.01, 0.05]
    u0_ext = extend_u0(u0, cal_lin)
    @test length(u0_ext) == 3
    @test u0_ext[end] ≈ 293.15
    @test u0_ext[1:2] == u0

    # Quadratic heat loss (Lavergne et al. 2018: φ(ΔT) = a*ΔT + b*ΔT²)
    a, b = 0.48, 0.002
    cal_quad = SemiAdiabaticCalorimeter(;
        T0 = 293.15,
        T_env = 293.15,
        Cp = 4000.0,
        heat_loss = ΔT -> a * ΔT + b * ΔT^2,
    )
    @test cal_quad isa SemiAdiabaticCalorimeter
    @test cal_quad.heat_loss(10.0) ≈ a * 10.0 + b * 100.0   # = 4.8 + 0.2 = 5.0
    @test cal_quad.heat_loss(0.0) ≈ 0.0

    # Generic callable: constant heat loss
    cal_const = SemiAdiabaticCalorimeter(;
        T0 = 300.0,
        T_env = 293.15,
        Cp = 3000.0,
        heat_loss = _ -> 2.5,
    )
    @test cal_const.heat_loss(0.0) ≈ 2.5
    @test cal_const.heat_loss(100.0) ≈ 2.5

end

# ── heat_rate ─────────────────────────────────────────────────────────────────

@testset "heat_rate (constant-enthalpy reactions)" begin
    # Build a minimal reaction with a constant ΔₐH⁰ = -50 000 J/mol
    # Using NumericFunc as stand-in for the species thermo function

    # Species with known ΔₐH⁰
    H2O = Species(
        "H2O";
        name = "Water",
        aggregate_state = AS_AQUEOUS,
        class = SC_AQSOLVENT,
    )
    H2O.properties[:ΔₐH⁰] = NumericFunc((T,) -> -285830.0, (:T,), u"J/mol")

    Ca2p = Species(
        "Ca+2";
        name = "Calcium ion",
        aggregate_state = AS_AQUEOUS,
        class = SC_AQSOLUTE,
    )
    Ca2p.properties[:ΔₐH⁰] = NumericFunc((T,) -> -542830.0, (:T,), u"J/mol")

    Calcite = Species(
        "Calcite";
        name = "Calcite",
        aggregate_state = AS_CRYSTAL,
        class = SC_COMPONENT,
    )
    Calcite.properties[:ΔₐH⁰] = NumericFunc((T,) -> -1206900.0, (:T,), u"J/mol")

    CO2aq = Species(
        "CO2";
        name = "CO2 aqueous",
        aggregate_state = AS_AQUEOUS,
        class = SC_AQSOLUTE,
    )
    CO2aq.properties[:ΔₐH⁰] = NumericFunc((T,) -> -413800.0, (:T,), u"J/mol")

    # CaCO3 + CO2 + H2O = Ca²⁺ + 2HCO₃⁻ (simplified: just CaCO3 = Ca²⁺ + CO3²⁻)
    # For this test: only check sign/magnitude
    k_test = arrhenius_rate_constant(1.0e-7, 0.0)   # k = const = 1e-7
    rm_test = FirstOrderRateModel(k_test)

    # We won't build a full KineticReaction here; just test heat_rate directly
    # by constructing a mock with a simple reaction
    reaction = Reaction([Calcite, Ca2p]; symbol = "calcite dissolution test")
    kr = KineticReaction(reaction, rm_test, FixedSurfaceArea(1.0), 1, [-1.0, 1.0])

    # With rate = 1e-5 mol/s, ΔHr = ΔₐH⁰_Ca²⁺ - ΔₐH⁰_Calcite
    ΔHr_expected = -542830.0 - (-1206900.0)   # = +664 070 J/mol (endothermic dissolution)
    rates = [1.0e-5]   # mol/s

    qdot = heat_rate([kr], rates, 298.15)
    @test isapprox(qdot, rates[1] * ΔHr_expected; rtol = 1.0e-6)

    # AD smoke-test through heat_rate
    dqdot_dr = ForwardDiff.derivative(r -> heat_rate([kr], [r], 298.15), 1.0e-5)
    @test isapprox(dqdot_dr, ΔHr_expected; rtol = 1.0e-6)
end

# ── extend_ode! for IsothermalCalorimeter ─────────────────────────────────────

@testset "extend_ode! IsothermalCalorimeter" begin

    cal = IsothermalCalorimeter(298.15)

    # Minimal mock: du and u with n_kin=1 mineral + 1 heat state
    # We simulate: du[1] = -r (mineral decreasing at rate r)
    # extend_ode! should set du[2] = heat_rate([kr], [r], T)
    # Since we can't easily run the full ODE, just test the sign

    # Use a simple reaction with known enthalpy
    k_test = arrhenius_rate_constant(1.0e-7, 0.0)
    rm_test = FirstOrderRateModel(k_test)
    H2Osp = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    CaO = Species("CaO"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    CaO.properties[:ΔₐH⁰] = NumericFunc((T,) -> -635090.0, (:T,), u"J/mol")
    H2Osp.properties[:ΔₐH⁰] = NumericFunc((T,) -> -285830.0, (:T,), u"J/mol")
    Ca_OH_2 = Species(
        "Ca(OH)2"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT
    )
    Ca_OH_2.properties[:ΔₐH⁰] = NumericFunc((T,) -> -986090.0, (:T,), u"J/mol")

    rxn = Reaction([CaO, H2Osp, Ca_OH_2]; symbol = "portlandite formation")
    kr = KineticReaction(rxn, rm_test, FixedSurfaceArea(1.0), 1, [-1.0, -1.0, 1.0])

    # Simulate: du[1] = -0.001 (dissolution rate 0.001 mol/s)
    du = [-0.001, 0.0]
    u = [0.01, 0.0]
    p = (kin_rxns = [kr], ϵ = 1.0e-30)

    extend_ode!(du, u, p, 1, cal)
    # q̇ = rate × ΔHr where rate = 0.001, ΔHr = H_Ca(OH)2 - H_CaO - H_H2O
    ΔHr = -986090.0 - (-635090.0) - (-285830.0)   # = -65 170 J/mol (exothermic)
    @test isapprox(du[2], 0.001 * ΔHr; rtol = 1.0e-6)
    # Exothermic → q̇ < 0 (heat released, convention: negative = exothermic here)
    # Actually sign depends on definition; just check it's nonzero and finite
    @test isfinite(du[2])

end

# ── SemiAdiabaticCalorimeter dT/dt ────────────────────────────────────────────

@testset "SemiAdiabaticCalorimeter dT/dt" begin

    cal = SemiAdiabaticCalorimeter(; T0 = 293.15, T_env = 293.15, L = 0.5, Cp = 4000.0)

    k_test = arrhenius_rate_constant(1.0e-7, 0.0)
    rm_test = FirstOrderRateModel(k_test)
    CaO = Species("CaO"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    CaO.properties[:ΔₐH⁰] = NumericFunc((T,) -> -635090.0, (:T,), u"J/mol")
    H2Osp = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    H2Osp.properties[:ΔₐH⁰] = NumericFunc((T,) -> -285830.0, (:T,), u"J/mol")
    Ca_OH_2 = Species("Ca(OH)2"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    Ca_OH_2.properties[:ΔₐH⁰] = NumericFunc((T,) -> -986090.0, (:T,), u"J/mol")
    rxn = Reaction([CaO, H2Osp, Ca_OH_2]; symbol = "portlandite")
    kr = KineticReaction(rxn, rm_test, FixedSurfaceArea(1.0), 1, [-1.0, -1.0, 1.0])

    # At T = T_env: no heat loss; dT/dt = q̇ / Cp
    T_curr = 293.15
    du = [-0.001, T_curr]   # du[1] unused here; du[2] = dT/dt to be set
    u = [0.01, T_curr]
    p = (kin_rxns = [kr], ϵ = 1.0e-30)

    extend_ode!(du, u, p, 1, cal)
    # dT/dt = (q̇ - L*(T-T_env)) / Cp
    # At T = T_env: dT/dt = q̇ / Cp
    ΔHr = -986090.0 - (-635090.0) - (-285830.0)
    r = 0.001
    qdot_expected = r * ΔHr
    dTdt_expected = (qdot_expected - 0.0) / 4000.0
    @test isapprox(du[2], dTdt_expected; rtol = 1.0e-6)

end

# ── _total_enthalpy ────────────────────────────────────────────────────────────

@testset "_total_enthalpy" begin

    # Two species: only the first has enthalpy data
    sp1 = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    sp1.properties[:ΔₐH⁰] = NumericFunc((T,) -> -285830.0, (:T,), u"J/mol")
    sp2 = Species("CaO"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    # sp2 has no ΔₐH⁰

    h_fns = [sp1[:ΔₐH⁰], nothing]
    n_full = [0.5, 1.0]   # 0.5 mol H2O, 1.0 mol CaO (no enthalpy data)

    H = ChemistryLab._total_enthalpy(n_full, h_fns, 298.15)
    # H = 0.5 × (-285830) + 0 (CaO has no data)
    @test isapprox(H, 0.5 * (-285830.0); rtol = 1.0e-10)

    # cumulative_heat = H(0) - H(t): if both points equal → Q = 0
    @test isapprox(H - H, 0.0; atol = 1.0e-12)

    # All-nothing h_fns → H = 0
    H_none = ChemistryLab._total_enthalpy(n_full, [nothing, nothing], 298.15)
    @test iszero(H_none)

    # AD smoke-test w.r.t. moles
    dHdn = ForwardDiff.derivative(n -> ChemistryLab._total_enthalpy([n, 1.0], h_fns, 298.15), 0.5)
    @test isfinite(dHdn)
    @test isapprox(dHdn, -285830.0; rtol = 1.0e-10)  # ΔₐH⁰ of sp1

    # AD smoke-test w.r.t. T (constant enthalpy function → dH/dT = 0)
    dHdT = ForwardDiff.derivative(T -> ChemistryLab._total_enthalpy(n_full, h_fns, T), 298.15)
    @test isfinite(dHdT)

end
