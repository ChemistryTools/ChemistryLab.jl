# The coupled run: `PoreHumidity` inside a rate law, integrated in time.
#
# 0.16.0 built `PoreHumidity` and wired `_humidity_at` to dispatch on it, and
# nothing ever integrated it. Doing so exposed a defect that no unit test on the
# object itself could reach: at full saturation a van Genuchten law has an
# unbounded `dp_c/dS`, so the differentiated rate law handed the implicit solver
# an infinite Jacobian entry and the integration did not advance past `t = 0`.
# These tests would have caught it, so they exist.

@testsection "a rate law reading the pore humidity integrates" begin
    subs = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    sel = split(
        "C3S C2S C3A C4AF Portlandite Jennite ettringite monosulphate12 " *
            "C3AH6 C3FH6 H2O@"
    )
    cs = ChemicalSystem(
        speciation(subs, sel; aggregate_state = [AS_AQUEOUS]), CEMDATA_PRIMARIES
    )
    sp(n) = cs[n]
    law = VanGenuchten(; a = 37.5479e6, m = 1 / 2.1684)     # Baroghel-Bouny mix CO

    function build(wc; humidity)
        st = ChemicalState(cs)
        for (k, v) in pairs((C3S = 0.619, C2S = 0.165, C3A = 0.08, C4AF = 0.087))
            set_quantity!(st, string(k), v * u"kg")
        end
        set_quantity!(st, "H2O@", wc * u"kg")
        h = humidity === :pore ? PoreHumidity(law, cs; reference = st) : humidity
        rxns = AbstractReaction[]
        specs = (
            (
                "C3S", PK84_PARAMS_C3S,
                OrderedDict(sp("C3S") => 1.0, sp("H2O@") => 103 / 30),
                OrderedDict(sp("Jennite") => 1.0, sp("Portlandite") => 4 / 3),
            ),
            (
                "C2S", PK84_PARAMS_C2S,
                OrderedDict(sp("C2S") => 1.0, sp("H2O@") => 73 / 30),
                OrderedDict(sp("Jennite") => 1.0, sp("Portlandite") => 1 / 3),
            ),
            (
                "C3A", PK84_PARAMS_C3A,
                OrderedDict(sp("C3A") => 1.0, sp("H2O@") => 6.0),
                OrderedDict(sp("C3AH6") => 1.0),
            ),
            (
                "C4AF", PK84_PARAMS_C4AF,
                OrderedDict(
                    sp("C4AF") => 1.0, sp("Portlandite") => 2.0, sp("H2O@") => 10.0
                ),
                OrderedDict(sp("C3AH6") => 1.0, sp("C3FH6") => 1.0),
            ),
        )
        for (nm, pk, r, p) in specs
            rx = Reaction(r, p; symbol = nm)
            rx[:rate] = parrot_killoh_avrami(
                pk, nm; α_max = 1.0, blaine = 380.0u"m^2/kg", humidity = h
            )
            push!(rxns, rx)
        end
        kp = KineticsProblem(cs, rxns, st, (0.0, 90 * 86400.0))
        ks = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-10)
        return kp, integrate(kp, ks), h
    end

    α_end(sol, kp) = let d = degrees_of_hydration(sol, kp; times = [sol.t[end]])
        v = d["C3S"]
        v isa AbstractVector ? v[end] : v
    end

    # ── the integration reaches the end of its span, which is the regression ──
    kp_free, sol_free, _ = build(0.35; humidity = nothing)
    kp_pore, sol_pore, h = build(0.35; humidity = :pore)
    @test sol_free.t[end] ≈ 90 * 86400.0
    @test sol_pore.t[end] ≈ 90 * 86400.0     # was 0.0 before the singular point
    # at S = 1 was regularized

    # ── the arrest is a result, not an imposed cap ────────────────────────────
    a_free, a_pore = α_end(sol_free, kp_free), α_end(sol_pore, kp_pore)
    @test 0 < a_pore < a_free                # self-desiccation stops it earlier
    @test a_free > 0.85                      # and without it the reaction runs on

    # ── it stops where `humidity_factor` cuts, and not somewhere else ─────────
    n_end = ustrip.(us"mol", state_at(sol_pore, kp_pore, sol_pore.t[end]).n)
    @test isapprox(h(n_end), 0.8; atol = 0.02)
    # The saturation at arrest is the one the static water budget uses: the
    # retention law puts RH 0.80 at S* = 0.786, and the trajectory lands there
    # on its own. Two independent routes to the same number.
    @test isapprox(pore_saturation(h, n_end), 0.786; atol = 0.01)

    # ── more water, more hydration ────────────────────────────────────────────
    kp_wet, sol_wet, _ = build(0.45; humidity = :pore)
    @test α_end(sol_wet, kp_wet) > a_pore

    # ── no regression for the two humidity forms that worked before ───────────
    kp_c, sol_c, _ = build(0.35; humidity = 0.95)          # a constant
    @test sol_c.t[end] ≈ 90 * 86400.0
    @test α_end(sol_c, kp_c) > 0
    kp_t, sol_t, _ = build(0.35; humidity = t -> 0.95)     # a function of time
    @test sol_t.t[end] ≈ 90 * 86400.0
    @test isapprox(α_end(sol_t, kp_t), α_end(sol_c, kp_c); rtol = 1.0e-6)
end

@testsection "the saturated point is regularized, not left infinite" begin
    # The unit-level statement of the same defect: the gradient of the humidity
    # with respect to the composition must be finite at full saturation, since
    # that is the initial condition of every sealed paste.
    subs = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    sel = split("C3S C2S Portlandite Jennite H2O@")
    cs = ChemicalSystem(
        speciation(subs, sel; aggregate_state = [AS_AQUEOUS]), CEMDATA_PRIMARIES
    )
    st = ChemicalState(cs)
    set_quantity!(st, "C3S", 0.8u"kg")
    set_quantity!(st, "C2S", 0.2u"kg")
    set_quantity!(st, "H2O@", 0.4u"kg")
    h = PoreHumidity(VanGenuchten(; a = 37.5479e6, m = 1 / 2.1684), cs; reference = st)
    n = ustrip.(us"mol", st.n)

    @test pore_saturation(h, n) ≈ 1.0
    @test h(n) ≈ 1.0
    g = ForwardDiff.gradient(x -> h(x), n)
    @test all(isfinite, g)                   # was an infinite norm
    @test all(iszero, g)                     # the regularization, stated

    # Just off saturation the real derivative is used, and it is not zero.
    n2 = copy(n)
    i_w = only(cs.idx_solvent)
    n2[i_w] *= 0.9
    @test h(n2) < 1.0
    g2 = ForwardDiff.gradient(x -> h(x), n2)
    @test all(isfinite, g2)
    @test any(!iszero, g2)
end
