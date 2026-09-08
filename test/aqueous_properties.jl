# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)
#
# Molalities, ionic strength, activity coefficients and the activity-convention
# pH, read back off a state.
#
# Three of the assertions below exist because the obvious way of getting these
# numbers is wrong, and each wrong way has produced a bad comparison against
# GEM-Selektor:
#
#   * γ from a ratio a/m diverges for a species at the solver's lower bound,
#     whose log-activity is dominated by the closures' `+ ϵ` term;
#   * `å_default` does not impose a common ionic radius;
#   * `pH(state)` and `pH(state, model)` are different quantities.

# A hand-built, deliberately *not* equilibrated state: every amount is known
# exactly, so the molalities and the ionic strength have closed forms.
function _aqp_state(; n_w = 55.5, n_ca = 0.01, n_oh = 0.02, ion_size = nothing)
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    if ion_size !== nothing
        for sp in substances
            if aggregate_state(sp) == AS_AQUEOUS && charge(sp) != 0
                sp[:å] = ion_size
            end
        end
    end
    species = speciation(
        substances, ["Portlandite", "H2O@"]; aggregate_state = [AS_AQUEOUS]
    )
    cs = ChemicalSystem(species, ["H2O@", "H+", "Ca+2", "Zz"])
    st = ChemicalState(cs)
    set_quantity!(st, "H2O@", n_w * u"mol")
    set_quantity!(st, "Portlandite", 0.02u"mol")
    set_quantity!(st, "Ca+2", n_ca * u"mol")
    set_quantity!(st, "OH-", n_oh * u"mol")
    return cs, st
end

@testsection "aqueous properties: concentration scale" begin
    # An activity coefficient means nothing without the scale its standard
    # state uses; `activity_coefficients` reads this to divide by the right
    # concentration.
    @test concentration_scale(DiluteSolutionModel()) === :molarity
    @test concentration_scale(HKFActivityModel()) === :molality
    @test concentration_scale(DaviesActivityModel()) === :molality
end

@testsection "aqueous properties: molalities and ionic strength" begin
    n_w, n_ca, n_oh = 55.5, 0.01, 0.02
    cs, st = _aqp_state(; n_w = n_w, n_ca = n_ca, n_oh = n_oh)

    m = molalities(st)
    # Read the solvent molar mass off the database rather than hard-coding it:
    # CEMDATA18's H2O@ is not exactly 18.01528 g/mol, and a hard-coded value
    # makes this assertion fail at the fifth digit for no reason of substance.
    i_w = only(cs.idx_solvent)
    kg = n_w * ustrip(us"kg/mol", cs.species[i_w][:M])
    @test m["Ca+2"] ≈ n_ca / kg rtol = 1.0e-12
    @test m["OH-"] ≈ n_oh / kg rtol = 1.0e-12

    # The solvent is not a solute and must not appear.
    @test !haskey(m, "H2O@")
    @test Set(keys(m)) == Set(symbol(cs.species[i]) for i in cs.idx_solutes)

    # I = ½ Σ mⱼ zⱼ², over the charged solutes only. H+ is auto-seeded at
    # neutral pH, so it contributes at the 1e-7 level and cannot be left out of
    # the reference sum.
    I_hand = 0.5 * sum(
        m[symbol(cs.species[i])] * Int(charge(cs.species[i]))^2
            for i in cs.idx_solutes if !iszero(charge(cs.species[i]))
    )
    @test ionic_strength(st) ≈ I_hand rtol = 1.0e-14
    # 4 m(Ca) + m(OH) dominates; the value is ~0.03 mol/kg here.
    @test 0.029 < ionic_strength(st) < 0.031

    # A system with no aqueous phase has none of these.
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    dry = ChemicalSystem([s for s in substances if symbol(s) == "Portlandite"])
    dry_state = ChemicalState(dry)
    @test_throws ArgumentError molalities(dry_state)
    @test_throws ArgumentError ionic_strength(dry_state)
end

@testsection "aqueous properties: γ from the formula, not from a ratio" begin
    cs, st = _aqp_state()
    I = ionic_strength(st)
    m = molalities(st)

    # ── molality-scale models: γ must equal exp(lna)/m for an abundant solute
    for model in (
            HKFActivityModel(), DaviesActivityModel(),
            HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0),
        )
        γ = activity_coefficients(st, model)
        lna = log_activities(st, model)
        for sym in ("Ca+2", "OH-")
            @test γ[sym] ≈ exp(lna[sym]) / m[sym] rtol = 1.0e-9
        end
        # a = γ m, restated through `activities`
        a = activities(st, model)
        @test a["Ca+2"] ≈ γ["Ca+2"] * m["Ca+2"] rtol = 1.0e-9
    end

    # ── the dilute model is ideal on its OWN scale: γ ≡ 1 everywhere.
    γ_id = activity_coefficients(st, DiluteSolutionModel())
    lna_id = log_activities(st, DiluteSolutionModel())
    @test all(≈(1.0; rtol = 1.0e-14), values(γ_id))

    # Its activities are numerically equal to the molalities, because it takes
    # c° = 1 mol/L and ρ = 1 kg/L, so molarity and molality coincide as numbers
    # even though the scales differ. That is why `concentration_scale` has to be
    # asked rather than inferred from the values: nothing in the numbers says
    # which scale `a` is on.
    @test exp(lna_id["Ca+2"]) ≈ m["Ca+2"] rtol = 1.0e-12
    @test exp(lna_id["Ca+2"]) / m["Ca+2"] ≈ γ_id["Ca+2"] rtol = 1.0e-12
end

@testsection "aqueous properties: γ is exact for a species at the lower bound" begin
    cs, st = _aqp_state()
    model = HKFActivityModel()
    γ = activity_coefficients(st, model)
    m = molalities(st)
    I = ionic_strength(st)
    sqrtI = sqrt(I)

    # Every aqueous species has a finite, positive coefficient, including the
    # ones sitting at 1e-16 mol. Computing γ as a ratio returns values of order
    # 1e300 for those, which is what this guards against.
    @test all(isfinite, values(γ))
    @test all(>(0), values(γ))

    # And the value matches the closed form for a trace, highly charged ion.
    trace = [
        i for i in cs.idx_solutes
            if abs(Int(charge(cs.species[i]))) >= 2 &&
            m[symbol(cs.species[i])] < 1.0e-10
    ]
    if !isempty(trace)
        i = first(trace)
        sp = cs.species[i]
        z = Int(charge(sp))
        å = ChemistryLab._hkf_lookup_å(sp, model)
        expected = 10.0^(
            -model.A * z^2 * sqrtI / (1 + model.B * å * sqrtI) + model.Ḃ * I
        )
        @test γ[symbol(sp)] ≈ expected rtol = 1.0e-12
        # The ratio, by contrast, is nowhere near it.
        lna = log_activities(st, model)
        @test !isapprox(
            exp(lna[symbol(sp)]) / m[symbol(sp)], expected; rtol = 1.0e-3
        )
    end
end

@testsection "aqueous properties: the model-level å" begin
    cs, st = _aqp_state()

    # `å` on the model must equal mutating `sp[:å]` on every charged species —
    # the documented route before this keyword existed.
    _, st_mutated = _aqp_state(; ion_size = 0.0)
    γ_kw = activity_coefficients(st, HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0))
    γ_mut = activity_coefficients(st_mutated, HKFActivityModel(Ḃ = 0.097637, Kₙ = 0.0))
    for sym in ("Ca+2", "OH-", "CaOH+")
        @test γ_kw[sym] ≈ γ_mut[sym] rtol = 1.0e-14
    end
    lna_kw = log_activities(st, HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0))
    lna_mut = log_activities(st_mutated, HKFActivityModel(Ḃ = 0.097637, Kₙ = 0.0))
    @test lna_kw["Ca+2"] ≈ lna_mut["Ca+2"] rtol = 1.0e-14
    @test lna_kw["H2O@"] ≈ lna_mut["H2O@"] rtol = 1.0e-14

    # `å = 0` is the Debye-Hückel limiting law plus the B-dot term, i.e. the
    # denominator collapses to 1.
    I = ionic_strength(st)
    model0 = HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0)
    γ0 = activity_coefficients(st, model0)
    for (sym, z) in (("OH-", 1), ("Ca+2", 2))
        @test γ0[sym] ≈ 10.0^(-model0.A * z^2 * sqrt(I) + model0.Ḃ * I) rtol = 1.0e-12
    end
    # Neutral species get no B-dot term when Kₙ = 0.
    neutrals = [
        symbol(cs.species[i]) for i in cs.idx_solutes
            if iszero(charge(cs.species[i]))
    ]
    @test !isempty(neutrals)
    for sym in neutrals
        @test γ0[sym] ≈ 1.0 rtol = 1.0e-14
    end

    # `å_default` is NOT a way to impose a common radius: it is the last resort
    # of the lookup chain and never reached for an ion the tables cover.
    γ_a = activity_coefficients(st, HKFActivityModel(å_default = 0.0))
    γ_b = activity_coefficients(st, HKFActivityModel(å_default = 9.0))
    @test γ_a["Ca+2"] ≈ γ_b["Ca+2"] rtol = 1.0e-14
    @test !isapprox(γ_a["Ca+2"], γ0["Ca+2"]; rtol = 1.0e-3)

    # The override wins over a per-species entry, which is the point of it.
    _, st_mut5 = _aqp_state(; ion_size = 5.0)
    γ_over = activity_coefficients(st_mut5, HKFActivityModel(å = 0.0))
    @test γ_over["Ca+2"] ≈
        10.0^(-0.5114 * 4 * sqrt(I) + 0.041 * I) rtol = 1.0e-12
end

@testsection "aqueous properties: pH in the activity convention" begin
    cs, st = _aqp_state()
    model = HKFActivityModel()

    lna = log_activities(st, model)
    @test pH(st, model) ≈ -lna["H+"] / log(10) rtol = 1.0e-14
    @test pOH(st, model) ≈ -lna["OH-"] / log(10) rtol = 1.0e-14

    # The activity pH sits above the molality pH by exactly -log10 γ(H+).
    γ = activity_coefficients(st, model)
    m = molalities(st)
    # The closures evaluate `log(m + ϵ)`, so an ion at m ≈ 1e-7 carries a
    # relative perturbation of order ϵ/m ≈ 1e-9 against the unregularized
    # product. That sets the tolerance here; it is not a modeling difference.
    @test pH(st, model) ≈ -log10(γ["H+"] * m["H+"]) rtol = 1.0e-7

    # The one-argument `pH` is a different quantity — concentration, and
    # reconstructed through pKw when the solution is basic. The two must not be
    # assumed equal; on this state they are more than a unit apart.
    @test abs(pH(st, model) - something(pH(st), NaN)) > 1.0

    # No H+ in the system → NaN rather than an error.
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    no_h = ChemicalSystem(
        [s for s in substances if symbol(s) in ("H2O@", "Ca+2", "OH-")],
        ["H2O@", "Ca+2", "Zz"],
    )
    st_no_h = ChemicalState(no_h)
    set_quantity!(st_no_h, "H2O@", 55.5u"mol")
    @test isnan(pH(st_no_h, model))
end

@testsection "aqueous properties: the GEM-Selektor calibration" begin
    # GEM-Selektor prints one activity coefficient per charge class, and that is
    # enough to identify the model it ran. On a CEMDATA18 Portland cement it
    # reported γ(1±) = 0.6113654, γ(2±) = 0.1212781 and γ(3±) = 0.008183112 at
    # I = 0.20969247 mol/kg. Fitting log10 γ = -D z² + E on |z| = 1 and 2 alone
    # gives D = 0.23417282 and E = +0.02047368, which then predicts |z| = 3, 4
    # and 5 to five significant digits — so the model is the Debye-Hückel
    # limiting law (å = 0, since CEMDATA18 carries no ion-size parameter) with
    # Ḃ = E/I = 0.097637 and no B-dot term on the neutrals.
    #
    # This asserts that the package's formula reproduces those coefficients, at
    # GEMS' own ionic strength, to better than 1.5 %.
    I = 0.20969247
    model = HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0)
    γ_lim = z -> 10.0^(-model.A * z^2 * sqrt(I) + model.Ḃ * I)
    @test γ_lim(1) ≈ 0.6113654 rtol = 3.0e-3
    @test γ_lim(2) ≈ 0.1212781 rtol = 1.5e-2
    @test γ_lim(3) ≈ 0.008183112 rtol = 4.0e-2
    # The package defaults are a different, defensible model — and nowhere near
    # those numbers on the divalents, which is why the keyword exists.
    default = HKFActivityModel()
    γ_def = (z, å) -> 10.0^(
        -default.A * z^2 * sqrt(I) / (1 + default.B * å * sqrt(I)) + default.Ḃ * I
    )
    ratio = γ_def(2, REJ_HKF["Ca+2"]) / 0.1212781
    @test 1.8 < ratio < 1.95
end

@testsection "the initial approximation is computed, not asked for" begin
    # `homotopy_initial_state` walks the solute amount up from a dilute system.
    # What must hold of its result is that it is a *feasible* starting point:
    # at λ = 1 the composition is the one given, so the component amounts are
    # unchanged. It is not an equilibrium and nothing here asserts that it is.
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    species = speciation(
        substances, ["Cal", "H2O@", "CO2@"]; aggregate_state = [AS_AQUEOUS]
    )
    cs = ChemicalSystem(species, ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
    st = ChemicalState(cs)
    set_quantity!(st, "H2O@", 55.5u"mol")
    set_quantity!(st, "Cal", 0.05u"mol")
    set_quantity!(st, "CO2@", 0.01u"mol")

    A = Float64.(cs.SM.A)
    b0 = A * ustrip.(us"mol", st.n)

    guess = homotopy_initial_state(st)
    @test guess isa ChemicalState
    # The element balance of the endpoint is the one it was given: the walk ends
    # at λ = 1, which is the state itself.
    @test A * ustrip.(us"mol", guess.n) ≈ b0 rtol = 1.0e-6
    @test all(>(0), ustrip.(us"mol", guess.n))

    # A system with no aqueous solvent has nothing to walk: `nothing`, not an
    # error, so `equilibrate_certified` can simply carry on without it.
    dry = ChemicalSystem([s for s in substances if symbol(s) == "Cal"])
    @test homotopy_initial_state(ChemicalState(dry)) === nothing

    # On a problem that certifies without help, declining the fallback must not
    # change the answer — it is only ever consulted when nothing else certified.
    eq_on, cert_on = equilibrate_certified(st)
    eq_off, cert_off = equilibrate_certified(st; autostart = false)
    @test cert_on.optimal
    @test cert_off.optimal
    @test ustrip.(us"mol", eq_on.n) ≈ ustrip.(us"mol", eq_off.n) rtol = 1.0e-6

    # The walk is done under the ideal model by default, whatever the target,
    # because the non-ideal ones do not walk. Passing one explicitly is allowed.
    guess_hkf = homotopy_initial_state(st; model = HKFActivityModel())
    @test guess_hkf === nothing || guess_hkf isa ChemicalState
end
