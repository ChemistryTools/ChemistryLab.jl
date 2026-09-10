# The capillary coupling: water retention, the Kelvin relation, the constraint
# that shifts the solvent's chemical potential, and the humidity a rate law can
# read off the state.
#
# The one thing this file must not do is suggest that the Kelvin term arrests
# hydration. It does not, by two orders of magnitude, and the last section pins
# that measurement so nobody later "improves" the package by claiming otherwise.

@testsection "Kelvin, and the two conventions" begin
    γ = 0.0728u"N/m"        # water against its vapor at 20 °C
    V_m = 1.807e-5u"m^3/mol"
    T = 298.15u"K"

    # A round trip, and the figure that anchors the whole feature: the
    # self-desiccation plateau of a sealed high-performance paste sits at the
    # gel-pore scale.
    r = kelvin_radius(0.8; γ = γ, V_m = V_m, T = T)
    @test 4.0e-9 < r < 5.5e-9
    @test kelvin_activity(r * u"m"; γ = γ, V_m = V_m, T = T) ≈ 0.8 rtol = 1.0e-10

    # Monotone, and bounded by one: a wider pore holds its water less tightly,
    # and a flat meniscus holds it not at all.
    a_small = kelvin_activity(2.0e-9u"m"; γ = γ, V_m = V_m, T = T)
    a_large = kelvin_activity(50.0e-9u"m"; γ = γ, V_m = V_m, T = T)
    @test a_small < a_large < 1.0

    # The two limits are refused rather than returned as Inf or a NaN.
    @test_throws ArgumentError kelvin_radius(1.0; γ = γ, V_m = V_m, T = T)
    @test_throws ArgumentError kelvin_radius(0.0; γ = γ, V_m = V_m, T = T)
    @test_throws ArgumentError kelvin_activity(0.0u"m"; γ = γ, V_m = V_m, T = T)

    # Every parameter is mandatory. The surface tension belongs to the liquid and
    # the temperature; the molar volume to the liquid. Neither is the package's
    # to assume, and Julia enforces that for free.
    @test_throws UndefKeywordError kelvin_activity(4.0e-9u"m"; γ = γ, V_m = V_m)
    @test_throws UndefKeywordError kelvin_radius(0.8; γ = γ, T = T)

end

@testsection "water retention laws" begin
    # ── no defaults, and the language is what enforces it ────────────────────
    @test_throws UndefKeywordError VanGenuchten(; a = 37.5479e6)
    @test_throws UndefKeywordError VanGenuchten(; m = 0.46)
    @test_throws UndefKeywordError TabulatedRetention(; S = [0.5, 1.0])

    # ── the published pair, and the b = 1/m convention ───────────────────────
    # Baroghel-Bouny et al. (1999), their Table 5, mix CO: an ordinary cement
    # paste at W/C = 0.34. They write the expression with `b = 1/m`, so a `b`
    # from that literature enters here as `m = 1/b`. Getting that inversion
    # wrong is silent and changes the exponent, which is why it is pinned.
    co = VanGenuchten(; a = 37.5479e6, m = 1 / 2.1684)
    @test co.m ≈ 0.46117 atol = 1.0e-6
    @test capillary_pressure(co, 1.0) == 0.0          # a flat meniscus holds nothing
    @test capillary_pressure(co, 0.5) > capillary_pressure(co, 0.9) > 0
    # Against the closed form, evaluated independently.
    S = 0.7
    @test capillary_pressure(co, S) ≈
        37.5479e6 * (S^(-2.1684) - 1)^(1 - 1 / 2.1684) rtol = 1.0e-12

    # Ranges are checked in an inner constructor rather than trusted.
    @test_throws ArgumentError VanGenuchten(; a = -1.0, m = 0.5)
    @test_throws ArgumentError VanGenuchten(; a = 1.0e6, m = 1.5)
    @test_throws ArgumentError VanGenuchten(; a = 1.0e6, m = 0.0)

    # ── a measured isotherm ─────────────────────────────────────────────────
    t = TabulatedRetention(; S = [0.3, 0.5, 0.75, 1.0], a_w = [0.44, 0.66, 0.85, 1.0])
    @test t(0.3) ≈ 0.44
    @test t(1.0) ≈ 1.0
    @test t(0.5) < t(0.6) < t(0.75)             # monotone between the knots
    # Linear in `ln a_w`, which is the variable the coupling reads.
    @test t(0.625) ≈ exp((log(0.66) + log(0.85)) / 2) rtol = 1.0e-12

    # A table that holds its water more tightly as the pore space fills is not a
    # retention curve, and interpolating it would be worse than refusing it.
    @test_throws ArgumentError TabulatedRetention(; S = [0.3, 0.5], a_w = [0.7, 0.5])
    @test_throws ArgumentError TabulatedRetention(; S = [0.5, 0.3], a_w = [0.5, 0.7])
    @test_throws ArgumentError TabulatedRetention(; S = [0.3, 0.5], a_w = [0.5, 1.5])
    @test_throws ArgumentError TabulatedRetention(; S = [0.3], a_w = [0.5])
    @test_throws ArgumentError TabulatedRetention(; S = [0.3, 0.5], a_w = [0.5])

    # Outside the measured range it clamps and says so once. An isotherm measured
    # between S = 0.3 and 1 says nothing at S = 0.05, and extrapolating a
    # logarithm there produces confident absurdity.
    @test_logs (:warn,) t(0.05)
    @test t(0.05) == 0.44                          # already warned, now silent
    @test t(1.2) == 1.0

    # An arbitrary function is a law too — the escape hatch.
    f = FunctionRetention(S -> 0.5 + 0.5 * S)
    @test f(0.0) ≈ 0.5
    @test f(1.0) ≈ 1.0

end

@testsection "the capillary constraint" begin
    sp = Dict(
        symbol(s) => s for s in build_species(
                datapath("slop98-inorganic-thermofun.json")
            )
    )
    cs = ChemicalSystem(
        [sp[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
        ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
    )
    fresh = ChemicalState(cs)
    set_quantity!(fresh, "H2O@", 1.0u"kg")
    set_quantity!(fresh, "Cal", 0.05u"mol")

    # ── the saturation the whole coupling rests on ───────────────────────────
    flat = FunctionRetention(_ -> 1.0)
    h = PoreHumidity(flat, cs; reference = fresh)
    n0 = ustrip.(us"mol", fresh.n)
    # Everything is liquid plus a little calcite, and the reference is the state
    # itself, so the pore space is what the liquid does not share with the solid.
    V_ref = ustrip(us"m^3", volume(fresh).total)
    V_liq = ustrip(us"m^3", volume(fresh).liquid)
    V_sol = ustrip(us"m^3", volume(fresh).solid)
    @test pore_saturation(h, n0) ≈ V_liq / (V_ref - V_sol) rtol = 1.0e-10
    @test 0 < pore_saturation(h, n0) <= 1

    # ── what the constraint refuses, and why ─────────────────────────────────
    # A system with no solvent has no water to hold.
    dry = ChemicalSystem([sp["Cal"]])
    @test_throws ArgumentError PoreHumidity(flat, dry; reference = ChemicalState(dry))
    # A reference smaller than the solids leaves no pore space to act on.
    @test_throws ArgumentError equilibrate_certified(
        fresh; constraint = CapillaryWater(flat, 1.0e-30u"m^3"),
    )
    # A law that claims the water is held even at full saturation is not a
    # retention curve: at S = 1 the meniscus is flat.
    @test_throws ArgumentError equilibrate_certified(
        fresh; constraint = CapillaryWater(FunctionRetention(_ -> 1.5); reference = fresh),
    )

    # ── the inert limit: the no-regression assertion that matters most ───────
    # A law that never holds the water must give back the unconstrained answer,
    # to the last digit, with a zero shift. If this drifts, the constraint is
    # changing answers it has no business changing.
    eq0, c0 = equilibrate_certified(fresh)
    @test c0.optimal
    qr = Ref(Float64[])
    eq1, c1 = equilibrate_certified(
        fresh; constraint = CapillaryWater(flat; reference = fresh), parameters = qr,
    )
    @test c1.optimal
    @test c1.param_residual < 1.0e-8
    @test only(qr[]) ≈ 0.0 atol = 1.0e-8
    @test ustrip.(us"mol", eq1.n) ≈ ustrip.(us"mol", eq0.n) rtol = 1.0e-6
    # And the 0.15.2 solvent guard must NOT fire: this is still a solution.
    @test solvent_fraction(eq1) > SOLVENT_FRACTION_FLOOR

    # ── a shift that is actually applied ────────────────────────────────────
    held = FunctionRetention(_ -> 0.9)
    qr2 = Ref(Float64[])
    eq2, c2 = equilibrate_certified(
        fresh; constraint = CapillaryWater(held; reference = fresh), parameters = qr2,
    )
    @test c2.optimal                                  # certified WITH the constraint
    @test c2.param_residual < 1.0e-8                  # and the closure holds
    @test only(qr2[]) ≈ log(0.9) rtol = 1.0e-6       # the shift is the one asked for
    # The shift lands on the water, which is the point of `hq` — but it is NOT
    # visible through `log_activities`, which evaluates the activity model and
    # knows nothing of the constraint. The activity the solver saw is the
    # model's plus the shift, and a caller who forgets that reads the
    # unshifted one and concludes nothing happened.
    lna_chem = log_activities(eq2, DiluteSolutionModel())["H2O@"]
    @test exp(lna_chem + only(qr2[])) ≈ 0.9 rtol = 1.0e-5
    @test exp(lna_chem + only(qr2[])) <
        exp(log_activities(eq0, DiluteSolutionModel())["H2O@"])
    # The unshifted value barely moves, which is the trap stated as a number.
    @test abs(lna_chem) < 1.0e-4

    # The old certificate could not see any of this. Audited as an unconstrained
    # problem, the same answer fails by about the size of the shift — which is
    # why `optimality_certificate` had to learn about constraints, and why this
    # is pinned as an anti-test.
    des = DualEquilibriumSolver(cs, DiluteSolutionModel())
    blind = optimality_certificate(des, eq2)
    @test !blind.optimal
    # Orders of magnitude worse than the constrained audit of the same
    # composition — stated relatively, since the absolute size depends on the
    # system and a fixed threshold would be a machine-dependent choice.
    @test blind.stationarity > 1.0e4 * max(c2.stationarity, 1.0e-16)

end

@testsection "the Kelvin term does not arrest hydration" begin
    # This section exists to stop a plausible-sounding claim from creeping back
    # in. Measured on a CEM I paste, imposing a water activity anywhere from
    # 0.95 down to 0.05 leaves the equilibrium assemblage unchanged; the arrest
    # of a real paste at 75-80 % RH is transport and nucleation, not
    # thermodynamics, and it belongs in `humidity_factor`.
    γ = 0.0728u"N/m"
    V_m = 1.807e-5u"m^3/mol"
    T = 298.15u"K"
    RT = 8.31446261815324 * 298.15

    # Alite going to a C-S-H plus portlandite consumes about 3.3 mol of water per
    # mole, so at a_w = 0.80 the whole capillary contribution to the affinity is
    # 3.3 * RT * ln(0.80) — under 2 kJ/mol, against a hydration Gibbs energy of
    # order -100 kJ/mol.
    capillary = 3.3 * RT * log(0.8)
    @test -2500 < capillary < -1500
    @test abs(capillary) < 0.05 * 100_000

    # The activity that WOULD null it, and the pore radius that activity implies.
    a_w_needed = exp(-100_000 / (3.3 * RT))
    @test a_w_needed < 1.0e-5
    # Smaller than a water molecule: the mechanism is not available.
    @test kelvin_radius(a_w_needed; γ = γ, V_m = V_m, T = T) < 1.0e-9

end

# ── the shift reaches the accessors ─────────────────────────────────────────

@testsection "the capillary shift can be read back off a state" begin
    # 0.16.0 shipped the shift in the solver's parameter block and nothing that
    # read a state knew about it, so `log_activities` reported the chemical
    # water activity after a solve posed at 0.90. Both accessors now take it.
    subs = build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false)
    d = Dict(symbol(s) => s for s in subs)
    cs = ChemicalSystem(
        [d[s] for s in split("H2O@ H+ OH- Na+ Cl-")], ["H2O@", "H+", "Na+", "Cl-", "Zz"]
    )
    st = ChemicalState(cs)
    set_quantity!(st, "H2O@", 1.0u"kg")
    set_quantity!(st, "Na+", 0.1u"mol")
    set_quantity!(st, "Cl-", 0.1u"mol")
    set_quantity!(st, "H+", 1.0e-7u"mol")
    set_quantity!(st, "OH-", 1.0e-7u"mol")
    model = HKFActivityModel()

    chem = water_activity(st, model)
    @test 0.99 < chem < 1.0                        # a dilute solution, barely lowered

    # The shift multiplies the activity, so it adds to the logarithm.
    shift = log(0.9)
    held = water_activity(st, model; kelvin_shift = shift)
    @test isapprox(held, chem * 0.9; rtol = 1.0e-12)
    @test isapprox(
        log_activities(st, model; kelvin_shift = shift)["H2O@"],
        log_activities(st, model)["H2O@"] + shift; atol = 1.0e-14,
    )

    # A zero shift changes nothing, which is what makes the keyword safe to add.
    @test log_activities(st, model; kelvin_shift = 0.0) == log_activities(st, model)
    # Only the solvent moves.
    a0 = log_activities(st, model)
    a1 = log_activities(st, model; kelvin_shift = shift)
    for k in keys(a0)
        k == "H2O@" && continue
        @test a0[k] == a1[k]
    end
end
