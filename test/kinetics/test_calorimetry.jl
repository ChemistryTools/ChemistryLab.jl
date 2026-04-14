using ForwardDiff
using OrderedCollections

# ── IsothermalCalorimeter ─────────────────────────────────────────────────────

@testset "IsothermalCalorimeter" begin

    cal = IsothermalCalorimeter(298.15)
    @test cal isa IsothermalCalorimeter
    @test cal.T ≈ 298.15us"K"

    @test n_extra_states(cal) == 1

    u0 = [0.01, 0.05]
    u0_ext = extend_u0(u0, cal)
    @test length(u0_ext) == 3
    @test u0_ext[end] == 0.0
    @test u0_ext[1:2] == u0

end

# ── SemiAdiabaticCalorimeter ──────────────────────────────────────────────────

@testset "SemiAdiabaticCalorimeter" begin

    # Linear heat loss via L keyword
    cal_lin = SemiAdiabaticCalorimeter(;
        Cp = 4000.0u"J/K", T_env = 293.15u"K", L = 0.5u"W/K", T0 = 293.15u"K",
    )
    @test cal_lin isa SemiAdiabaticCalorimeter
    @test cal_lin.Cp ≈ 4000.0us"J/K"
    @test cal_lin.T_env ≈ 293.15us"K"
    @test cal_lin.T0 ≈ 293.15us"K"
    @test cal_lin.heat_loss(1.0) ≈ 0.5
    @test cal_lin.heat_loss(2.0) ≈ 1.0

    @test n_extra_states(cal_lin) == 1

    u0 = [0.01, 0.05]
    u0_ext = extend_u0(u0, cal_lin)
    @test length(u0_ext) == 3
    @test u0_ext[end] ≈ 293.15
    @test u0_ext[1:2] == u0

    # Quadratic heat loss (Lavergne et al. 2018)
    a, b = 0.48, 0.002
    cal_quad = SemiAdiabaticCalorimeter(;
        Cp = 4000.0u"J/K",
        T_env = 293.15u"K",
        heat_loss = ΔT -> a * ΔT + b * ΔT^2,
        T0 = 293.15u"K",
    )
    @test cal_quad isa SemiAdiabaticCalorimeter
    @test cal_quad.heat_loss(10.0) ≈ a * 10.0 + b * 100.0
    @test cal_quad.heat_loss(0.0) ≈ 0.0

    # Plain Real (SI) inputs
    cal_si = SemiAdiabaticCalorimeter(; Cp = 3500.0, T_env = 295.0, L = 0.3, T0 = 295.0)
    @test cal_si.Cp ≈ 3500.0us"J/K"
    @test cal_si.T_env ≈ 295.0us"K"

    # Constant heat loss
    cal_const = SemiAdiabaticCalorimeter(;
        Cp = 3000.0u"J/K",
        T_env = 293.15u"K",
        heat_loss = _ -> 2.5,
        T0 = 300.0u"K",
    )
    @test cal_const.T0 ≈ 300.0us"K"
    @test cal_const.heat_loss(0.0) ≈ 2.5
    @test cal_const.heat_loss(100.0) ≈ 2.5

    # Requires either heat_loss or L
    @test_throws ArgumentError SemiAdiabaticCalorimeter(;
        Cp = 4000.0, T_env = 293.15, T0 = 293.15,
    )

end

# ── heat_rate ─────────────────────────────────────────────────────────────────

@testset "heat_rate (constant-enthalpy reactions)" begin

    # Species with known ΔₐH⁰
    H2O = Species("H2O"; name = "Water", aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    H2O.properties[:ΔₐH⁰] = NumericFunc((T,) -> -285830.0, (:T,), u"J/mol")

    Ca2p = Species("Ca+2"; name = "Calcium ion", aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    Ca2p.properties[:ΔₐH⁰] = NumericFunc((T,) -> -542830.0, (:T,), u"J/mol")

    Calcite = Species("Calcite"; name = "Calcite", aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    Calcite.properties[:ΔₐH⁰] = NumericFunc((T,) -> -1206900.0, (:T,), u"J/mol")

    # Explicit reactants/products so complete_thermo_functions! gives correct ΔᵣH⁰
    reaction = Reaction(
        OrderedDict(Calcite => 1.0),
        OrderedDict(Ca2p => 1.0);
        symbol = "calcite dissolution test",
    )
    dummy_fn = KineticFunc((T, P, t, n, lna, n0) -> 0.0, NamedTuple(), u"mol/s")
    kr = KineticReaction(reaction, dummy_fn, 1, [-1.0, 1.0])

    # Thermodynamic ΔᵣH⁰ = ΔₐH⁰(Ca²⁺) − ΔₐH⁰(Calcite) = +664 070 J/mol (endothermic)
    # heat_rate uses −ΔᵣH⁰: negative for endothermic (heat absorbed from calorimeter)
    ΔHr_thermo = -542830.0 - (-1206900.0)   # = +664 070 J/mol
    rates = [1.0e-5]

    qdot = heat_rate([kr], rates, 298.15)
    @test isapprox(qdot, rates[1] * (-ΔHr_thermo); rtol = 1.0e-6)

    # With explicit heat_per_mol: should override stoichiometric path
    # heat_per_mol > 0 means heat generated (exothermic convention)
    kr_explicit = KineticReaction(reaction, dummy_fn, 1, [-1.0, 1.0]; heat_per_mol = 50_000.0)
    qdot_explicit = heat_rate([kr_explicit], rates, 298.15)
    @test isapprox(qdot_explicit, 1.0e-5 * 50_000.0; rtol = 1.0e-10)

    # AD smoke-test through heat_rate
    dqdot_dr = ForwardDiff.derivative(r -> heat_rate([kr], [r], 298.15), 1.0e-5)
    @test isapprox(dqdot_dr, -ΔHr_thermo; rtol = 1.0e-6)

end

# ── extend_ode! for IsothermalCalorimeter ─────────────────────────────────────

@testset "extend_ode! IsothermalCalorimeter" begin

    cal = IsothermalCalorimeter(298.15)

    CaO = Species("CaO"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    CaO.properties[:ΔₐH⁰] = NumericFunc((T,) -> -635090.0, (:T,), u"J/mol")
    H2Osp = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    H2Osp.properties[:ΔₐH⁰] = NumericFunc((T,) -> -285830.0, (:T,), u"J/mol")
    Ca_OH_2 = Species("Ca(OH)2"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    Ca_OH_2.properties[:ΔₐH⁰] = NumericFunc((T,) -> -986090.0, (:T,), u"J/mol")

    dummy_fn = KineticFunc((T, P, t, n, lna, n0) -> 0.0, NamedTuple(), u"mol/s")
    rxn = Reaction(
        OrderedDict(CaO => 1.0, H2Osp => 1.0),
        OrderedDict(Ca_OH_2 => 1.0);
        symbol = "portlandite formation",
    )
    kr = KineticReaction(rxn, dummy_fn, 1, [-1.0, -1.0, 1.0])

    du = [-0.001, 0.0]
    u = [0.01, 0.0]
    p = (kin_rxns = [kr], ϵ = 1.0e-30, rates_buf = [0.001])

    extend_ode!(du, u, p, 1, cal)
    # ΔᵣH⁰ = ΔₐH⁰(Ca(OH)₂) − ΔₐH⁰(CaO) − ΔₐH⁰(H₂O) = −65 170 J/mol (exothermic)
    # heat_rate uses −ΔᵣH⁰ > 0: positive qdot for exothermic (heat generated)
    ΔHr_thermo = -986090.0 - (-635090.0) - (-285830.0)   # = −65 170 J/mol
    @test isapprox(du[2], -0.001 * ΔHr_thermo; rtol = 1.0e-6)
    @test isfinite(du[2])

end

# ── SemiAdiabaticCalorimeter dT/dt ────────────────────────────────────────────

@testset "SemiAdiabaticCalorimeter dT/dt" begin

    cal = SemiAdiabaticCalorimeter(;
        Cp = 4000.0u"J/K", T_env = 293.15u"K", L = 0.5u"W/K", T0 = 293.15u"K",
    )

    CaO = Species("CaO"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    CaO.properties[:ΔₐH⁰] = NumericFunc((T,) -> -635090.0, (:T,), u"J/mol")
    H2Osp = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    H2Osp.properties[:ΔₐH⁰] = NumericFunc((T,) -> -285830.0, (:T,), u"J/mol")
    Ca_OH_2 = Species("Ca(OH)2"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    Ca_OH_2.properties[:ΔₐH⁰] = NumericFunc((T,) -> -986090.0, (:T,), u"J/mol")

    dummy_fn = KineticFunc((T, P, t, n, lna, n0) -> 0.0, NamedTuple(), u"mol/s")
    rxn = Reaction(
        OrderedDict(CaO => 1.0, H2Osp => 1.0),
        OrderedDict(Ca_OH_2 => 1.0);
        symbol = "portlandite",
    )
    kr = KineticReaction(rxn, dummy_fn, 1, [-1.0, -1.0, 1.0])

    # ΔᵣH⁰ = ΔₐH⁰(Ca(OH)₂) − ΔₐH⁰(CaO) − ΔₐH⁰(H₂O) = −65 170 J/mol (exothermic)
    # heat_rate = r × (−ΔᵣH⁰) > 0: positive for exothermic → T rises
    ΔHr_thermo = -986090.0 - (-635090.0) - (-285830.0)   # = −65 170 J/mol
    r = 0.001
    qdot_expected = r * (-ΔHr_thermo)   # = +65.17 W (heat generated)

    # At T = T_env: no heat loss (L * ΔT = 0); dT/dt = q̇ / Cp_total
    T_curr = 293.15
    du = [-0.001, T_curr]
    u = [0.01, T_curr]
    # p needs: kin_rxns, rates_buf, n_full, cp_fns
    p = (
        kin_rxns = [kr], ϵ = 1.0e-30, rates_buf = [0.001],
        n_full = [0.01], cp_fns = [nothing],
    )

    extend_ode!(du, u, p, 1, cal)
    # Cp_total = cal.Cp (all cp_fns are nothing)
    dTdt_expected = (qdot_expected - 0.0) / 4000.0   # positive → T rises
    @test isapprox(du[2], dTdt_expected; rtol = 1.0e-6)

    # With nonzero ΔT: heat loss L * ΔT reduces dT/dt
    T_hot = 303.15   # +10 °C above T_env
    du_hot = [-0.001, T_hot]
    u_hot = [0.01, T_hot]
    p_hot = (
        kin_rxns = [kr], ϵ = 1.0e-30, rates_buf = [0.001],
        n_full = [0.01], cp_fns = [nothing],
    )
    extend_ode!(du_hot, u_hot, p_hot, 1, cal)
    ΔT = T_hot - 293.15
    dTdt_hot = (qdot_expected - 0.5 * ΔT) / 4000.0
    @test isapprox(du_hot[2], dTdt_hot; rtol = 1.0e-6)

    # Variable Cp_total: add a Cp° function for the mineral
    cp_fn = NumericFunc((T,) -> 100.0, (:T,), u"J/(mol*K)")   # 100 J/(mol·K)
    p_cp = (
        kin_rxns = [kr], ϵ = 1.0e-30, rates_buf = [0.001],
        n_full = [0.01], cp_fns = [cp_fn],
    )
    du_cp = [-0.001, T_curr]
    u_cp = [0.01, T_curr]
    extend_ode!(du_cp, u_cp, p_cp, 1, cal)
    Cp_total_expected = 4000.0 + 0.01 * 100.0   # = 4001 J/K
    dTdt_cp = (qdot_expected - 0.0) / Cp_total_expected
    @test isapprox(du_cp[2], dTdt_cp; rtol = 1.0e-6)

end

# ── _total_enthalpy ────────────────────────────────────────────────────────────

@testset "_total_enthalpy" begin

    sp1 = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    sp1.properties[:ΔₐH⁰] = NumericFunc((T,) -> -285830.0, (:T,), u"J/mol")
    sp2 = Species("CaO"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)

    h_fns = [sp1[:ΔₐH⁰], nothing]
    n_full = [0.5, 1.0]

    H = ChemistryLab._total_enthalpy(n_full, h_fns, 298.15)
    @test isapprox(H, 0.5 * (-285830.0); rtol = 1.0e-10)
    @test isapprox(H - H, 0.0; atol = 1.0e-12)

    H_none = ChemistryLab._total_enthalpy(n_full, [nothing, nothing], 298.15)
    @test iszero(H_none)

    dHdn = ForwardDiff.derivative(n -> ChemistryLab._total_enthalpy([n, 1.0], h_fns, 298.15), 0.5)
    @test isfinite(dHdn)
    @test isapprox(dHdn, -285830.0; rtol = 1.0e-10)

    dHdT = ForwardDiff.derivative(T -> ChemistryLab._total_enthalpy(n_full, h_fns, T), 298.15)
    @test isfinite(dHdT)

end
