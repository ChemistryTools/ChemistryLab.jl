# The certified equilibrium path.
#
# Every assertion here is on a measured property of the solvers, not on a target
# value: that the certificate decides, that it is not fooled by a start which is
# itself infeasible, and that the multi-start route certifies cases no single
# back end does.

@testsection "Certified equilibrium" begin

    sp = Dict(
        symbol(s) => s for s in build_species(
                datapath("slop98-inorganic-thermofun.json")
            )
    )
    species = [sp[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")]
    cs = ChemicalSystem(species, ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
    A = Matrix{Float64}(cs.SM.A)

    function calcite(; θ = 25.0, nco2 = 0.0)
        st = ChemicalState(cs)
        set_quantity!(st, "Cal", 1.0e-3u"mol")
        set_quantity!(st, "H2O@", 1.0u"kg")
        nco2 > 0 && set_quantity!(st, "CO2@", nco2 * u"mol")
        V = volume(st)
        set_quantity!(st, "H+", 1.0e-4u"mol/L" * V.liquid)
        set_quantity!(st, "OH-", 1.0e-10u"mol/L" * V.liquid)
        set_temperature!(st, (273.15 + θ) * u"K")
        return st
    end

    # Element balance, each row against its own budget. An absolute ∞-norm hides
    # the small rows: at ‖An − b‖∞ = 3e-6 the water row (111 mol) is satisfied to
    # 1.8e-8 while the charge row (budget 1e-4) is out by 3 %.
    function balance_rel(st0, eq)
        b = A * [ustrip(us"mol", x) for x in st0.n]
        r = A * [ustrip(us"mol", x) for x in eq.n] - b
        keep = abs.(b) .> 1.0e-8
        return any(keep) ? maximum(abs.(r[keep] ./ b[keep])) : maximum(abs, r)
    end

    @testsection "the certificate decides, and the answer is reproducible" begin
        st = calcite()
        eq1, c1 = equilibrate_certified(calcite())
        eq2, c2 = equilibrate_certified(calcite())
        @test c1.optimal
        @test pH(eq1) == pH(eq2)          # bit-for-bit, no path dependence
        @test balance_rel(st, eq1) < 1.0e-10
        @test c1.worst_supersaturation < 0 || c1.worst_supersaturation == -Inf
    end

    @testsection "the interior point alone is wrong on this case" begin
        # Not a target value: the point is that the certified route is two orders
        # of magnitude better on a balance the old default reported as fine.
        st = calcite()
        eq_bar = equilibrate(calcite(), OptimaOptimizer())
        eq_cert, cert = equilibrate_certified(calcite())
        @test balance_rel(st, eq_bar) > 1.0e-3      # measured at 3.0e-2
        @test balance_rel(st, eq_cert) < 1.0e-8
        @test cert.optimal
        # And they disagree on the answer, not merely on the residual.
        @test abs(pH(eq_bar) - pH(eq_cert)) > 1.0e-3
    end

    @testsection "b is fixed once, not taken from each start" begin
        # A start that violates the balance shifts the component totals by its own
        # infeasibility. Recomputing `b` from each start therefore poses a
        # different problem per start, and the dual solve certifies the answer to
        # the shifted one: measured, two "certified" compositions 0.2 % apart on
        # dissolved calcium, which a convex problem with one minimum cannot have.
        st = calcite()
        b0 = A * [ustrip(us"mol", x) for x in st.n]
        des = DualEquilibriumSolver(cs, DiluteSolutionModel())
        bar = equilibrate(calcite(), OptimaOptimizer())
        from_raw = SciMLBase.solve(des, calcite(); b = b0)
        from_bar = SciMLBase.solve(des, bar; b = b0)
        # With the same b both routes must land on the same minimum.
        @test isapprox(pH(from_raw), pH(from_bar); atol = 1.0e-6)
        # Without it, the second solves a different problem — the shift being
        # exactly the interior point's own infeasibility.
        shifted = SciMLBase.solve(des, bar)
        @test abs(pH(shifted) - pH(from_raw)) > 1.0e-3
    end

    @testsection "equilibrate certifies by default, and can be told not to" begin
        st = calcite()
        @test balance_rel(st, equilibrate(calcite())) < 1.0e-8
        @test balance_rel(st, equilibrate(calcite(); certify = false)) > 1.0e-3
    end

    @testsection "the dual route needs an aqueous phase and H2O@" begin
        # It parameterizes the interior variables by the solvent's potential, so
        # both conditions are structural. `equilibrate_certified` returns
        # `nothing` for the certificate on such a system rather than claiming a
        # proof it cannot give.
        @test ChemistryLab._dual_applicable(cs)

        solid_only = ChemicalSystem(
            [Species("NaCl"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)],
            ["NaCl"],
        )
        @test !ChemistryLab._dual_applicable(solid_only)

        # Aqueous, but the solvent is not called `H2O@`.
        no_solvent = ChemicalSystem([sp["Ca+2"], sp["CO3-2"]], ["Ca+2", "CO3-2"])
        @test !ChemistryLab._dual_applicable(no_solvent)
    end


    @testsection "STRICT_CONVERGENCE is honored by the certified route" begin
        # The flag existed only on the interior-point retcode, so a caller who set
        # it — asking that a non-converged solve never pass as a result — still got
        # an uncertified answer back with a `@warn`. That answer can violate the
        # element balance by moles and still look like an ordinary `ChemicalState`:
        # measured on a cement paste, a balance off by 6.7 mol, every hydrate at
        # zero, and a table of amounts that reads as a result.
        #
        # `-b` is infeasible by construction: every component budget is negative and
        # every stoichiometric coefficient is non-negative, so no composition with
        # `n >= 0` can meet it and no route can certify.
        st = calcite()
        b_bad = -(A * [ustrip(us"mol", x) for x in st.n])

        strict = ChemistryLab.STRICT_CONVERGENCE[]
        try
            ChemistryLab.STRICT_CONVERGENCE[] = true
            @test_throws "no route produced a certifiable equilibrium" equilibrate_certified(
                st; b = b_bad, autostart = false
            )

            # The default stays a warning: the answer is still the best one found,
            # and `optimality_certificate` is there to audit it.
            ChemistryLab.STRICT_CONVERGENCE[] = false
            eq, cert = equilibrate_certified(st; b = b_bad, autostart = false)
            @test cert.optimal == false
            @test eq isa ChemicalState
        finally
            ChemistryLab.STRICT_CONVERGENCE[] = strict
        end

        # And the flag is left exactly as it was found.
        @test ChemistryLab.STRICT_CONVERGENCE[] == strict

    end
end

@testsection "ForwardDiff through the certified route" begin

    # The certified route solves in real arithmetic — an active set has a discrete
    # component, so the map `b ↦ n*(b)` is smooth only piecewise and pushing duals
    # through the search would be wrong. The derivative is attached at the answer,
    # from the optimality conditions.
    #
    # This has to be tested rather than assumed: the certified path converts the
    # component totals to `Float64` to fix them once, and before the dual branch
    # existed that silently dropped every partial.
    sp = Dict(
        symbol(s) => s for s in build_species(
                datapath("slop98-inorganic-thermofun.json")
            )
    )
    cs = ChemicalSystem(
        [sp[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
        ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
    )
    idx = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))

    function pH_of(x)
        n = Any[fill(0.0 * oneunit(x) * u"mol", length(cs.species))...]
        n[idx["H2O@"]] = 55.5 * oneunit(x) * u"mol"
        n[idx["Cal"]] = 0.05 * oneunit(x) * u"mol"
        n[idx["CO2@"]] = x * u"mol"
        return pH(equilibrate(ChemicalState(cs, n)))
    end

    x0 = 0.01
    ad = ForwardDiff.derivative(pH_of, x0)
    h = 1.0e-6
    fd = (pH_of(x0 + h) - pH_of(x0 - h)) / (2h)

    @test isfinite(ad)
    @test ad != 0                       # the partials are not dropped
    @test ad ≈ fd rtol = 1.0e-5         # measured at 4.3e-9
    @test ad < 0                        # more CO2, lower pH

    # And the primal value is the certified one, not something the dual path
    # computed separately.
    @test pH_of(x0) ≈ pH(
        first(
            equilibrate_certified(
                let
                    n = Any[fill(0.0u"mol", length(cs.species))...]
                    n[idx["H2O@"]] = 55.5u"mol"; n[idx["Cal"]] = 0.05u"mol"
                    n[idx["CO2@"]] = x0 * u"mol"
                    ChemicalState(cs, n)
                end
            )
        )
    ) atol = 1.0e-12

end

@testsection "the certified route restarts from its own answer" begin
    # `_keep_better` is the arbiter of every round of the multi-start search:
    # certified beats uncertified, and among uncertified the smaller KKT error
    # wins. It is unit-testable without a solve, and it is what keeps a restart
    # from ever making the answer worse.
    kb = ChemistryLab._keep_better
    a, b = :A, :B
    cert(opt, stat; bal = 0.0, si = -1.0) =
        (;
        optimal = opt, stationarity = stat, balance = bal,
        worst_supersaturation = si,
    )

    # Certified beats uncertified, in both directions and whatever the errors.
    @test first(kb(a, cert(false, 1.0e-16), b, cert(true, 1.0e-3))) === b
    @test first(kb(a, cert(true, 1.0e-3), b, cert(false, 1.0e-16))) === a

    # Among uncertified, the smaller KKT error wins.
    @test first(kb(a, cert(false, 1.0e-6), b, cert(false, 1.0e-9))) === b
    @test first(kb(a, cert(false, 1.0e-9), b, cert(false, 1.0e-6))) === a

    # And that error is the worst of ALL THREE residuals, not the stationarity
    # alone. This is the case that was wrong, and it is not academic: a Windows
    # run of a CEM I paste came back stationary to 2.4e-3 with an element
    # balance off by 6.7 mol and a phase supersaturated by 45, and it beat every
    # candidate the continuation produced — those being stationary to only 1e-2
    # while conserving mass. Ranked on stationarity, the answer that is not an
    # answer wins; ranked on the worst residual, the usable one does. It is also
    # the ranking `solve_certified` already used internally, so the two agree.
    bad = cert(false, 2.4e-3; bal = 6.7, si = 45.3)
    good = cert(false, 1.0e-2; bal = 1.0e-13, si = -0.5)
    @test ChemistryLab._kkt_error(bad) > ChemistryLab._kkt_error(good)
    @test first(kb(a, bad, b, good)) === b
    @test first(kb(a, good, b, bad)) === a

    # A negative worst supersaturation is not an error: every absent phase
    # undersaturated is what optimality requires, so it must not be counted.
    @test ChemistryLab._kkt_error(cert(false, 1.0e-9; si = -12.0)) == 1.0e-9

    # Between two certified answers the KKT error still decides, and a tie keeps
    # the incumbent: a round that buys nothing changes nothing, which is what
    # makes the restart loop safe to run and what stops it.
    @test first(kb(a, cert(true, 1.0e-12), b, cert(true, 1.0e-16))) === b
    @test first(kb(a, cert(false, 1.0e-9), b, cert(false, 1.0e-9))) === a
    @test first(kb(a, cert(true, 1.0e-16), b, cert(true, 1.0e-16))) === a

    # The bound exists so a case improving by a hair every round cannot loop.
    @test ChemistryLab._MAX_RESTARTS isa Integer
    @test ChemistryLab._MAX_RESTARTS >= 1

end

@testsection "a candidate's diagnostics are not the answer's" begin
    # `equilibrate_certified` runs every back end from several compositions
    # precisely because none of them works on every problem, and keeps whichever
    # answer the certificate proves. A candidate that does not converge is
    # therefore ordinary. Left unguarded it printed "returned `MaxIters`" and
    # "did not certify optimality" from candidates along the way, so a call that
    # ended `optimal = true` read as a failed solve.

    # The scope sets the flag and restores it, including when the body throws —
    # a start search that raises must not leave the whole session quiet.
    @test ChemistryLab._EXPLORING_STARTS[] == false
    inside = ChemistryLab._exploring_starts() do
        ChemistryLab._EXPLORING_STARTS[]
    end
    @test inside
    @test ChemistryLab._EXPLORING_STARTS[] == false
    @test_throws ErrorException ChemistryLab._exploring_starts() do
        error("a back end failed")
    end
    @test ChemistryLab._EXPLORING_STARTS[] == false

    # Nested scopes restore the previous value, not `false`: the walk runs inside
    # the search, and the inner scope ending must not un-quiet the outer one.
    ChemistryLab._exploring_starts() do
        ChemistryLab._exploring_starts() do
        end
        @test ChemistryLab._EXPLORING_STARTS[]
    end
    @test ChemistryLab._EXPLORING_STARTS[] == false

    # The contract as a caller sees it: a call that ends with a certificate emits
    # nothing at all. `@test_logs` installs a fresh logger, so the `maxlog = 1`
    # carried by those warnings does not make this depend on what ran before.
    #
    # This one assertion also guards the other silence in this release: SciMLBase
    # warns "arrays or dicts to store parameters of different types can hurt
    # performance" the moment a conservation matrix with an abstract element type
    # reaches the problem's parameters, so a regression in
    # `_concrete_conservation` shows up here as a stray record.
    # Built here rather than with `calcite()`, which is a local of another
    # testset in this file.
    sp3 = Dict(
        symbol(s) => s for s in build_species(
                datapath("slop98-inorganic-thermofun.json")
            )
    )
    cs3 = ChemicalSystem(
        [sp3[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
        ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
    )
    st3 = ChemicalState(cs3)
    set_quantity!(st3, "H2O@", 1.0u"kg")
    set_quantity!(st3, "Cal", 1.0e-3u"mol")
    eq3, cert3 = @test_logs equilibrate_certified(st3)
    @test cert3.optimal

end

@testsection "the conservation matrix reaching the solver is concretely typed" begin
    # `ChemicalSystem` stores its stoichiometry as `Matrix{Real}` whenever
    # integer and rational coefficients coexist — a cement's does, through
    # `C3AFS0.84H4.32` and its kind — which is right for the chemistry and wrong
    # for the solver. An abstract element type boxes every entry and makes
    # `mul!(res, A, x)` the generic fallback with a dispatch per element, on a
    # product evaluated at every objective and constraint call. SciMLBase warns
    # about exactly this as soon as such an array reaches a problem's parameters,
    # and the warning was correct.
    abstract_A = Matrix{Real}([1 0 2 // 5; 0 1 3])
    @test !isconcretetype(eltype(abstract_A))
    narrowed = ChemistryLab._concrete_float(abstract_A)
    @test isconcretetype(eltype(narrowed))
    @test eltype(narrowed) <: AbstractFloat      # not Rational: the pipeline is float
    @test !ChemistryLab.SciMLBase.should_warn_paramtype(narrowed)

    # Lossy for a non-dyadic rational, and deliberately so: `2//5` becomes `0.4`,
    # which is not equal to it. It is the same rounding the rest of the solver
    # already applies, and the exact matrix stays in `system.SM.A`.
    @test narrowed ≈ abstract_A
    @test narrowed[1, 3] == 0.4
    @test narrowed[1, 3] != 2 // 5

    # A concrete array is returned untouched, identically — the narrowing is for
    # the abstract case and nothing else. An exact rational stoichiometry stays
    # exact, and a caller differentiating through `A` keeps their number type.
    for m in (Rational{Int}[1 0; 0 1], Int[1 2; 3 4], Float32[1 0; 0 1])
        @test ChemistryLab._concrete_float(m) === m
    end

    # And what the solver actually receives. The default `b = A * u0` inherits
    # the abstract element type from `A`, so it has to be narrowed too — this is
    # the assertion that caught it.
    ep = EquilibriumProblem(abstract_A, (n, q) -> n, [1.0, 1.0, 1.0])
    @test isconcretetype(eltype(ep.A))
    @test isconcretetype(eltype(ep.b))
    @test ep.A ≈ abstract_A
    @test !ChemistryLab.SciMLBase.should_warn_paramtype(
        (; A = ep.A, b = ep.b, T = 298.15)
    )

end

@testsection "the certificate names the missing phase, and the route puts it in" begin
    # A positive worst supersaturation means a phase sits at the lower bound
    # while the solution is supersaturated with respect to it. On a convex
    # problem that is a genuine KKT failure: the active set is wrong and the
    # answer is not the answer. The certificate reports it as one number;
    # `saturation_indices` says which phase, and `_repair_start` puts it in.
    #
    # The case this exists for is a phase SWAP, which an active-set loop that
    # admits one phase at a time cannot perform. Reproduced exactly on a CEM I
    # paste by solving it with `hydrotalcite` out of the phase list: all
    # 0.02515 mol of magnesium goes to brucite, the aqueous phase equilibrates
    # with that assemblage, and putting hydrotalcite back leaves it absent and
    # supersaturated by 5.58 log units while the solve is otherwise impeccable —
    # stationarity 5.9e-16, element balance 1.6e-14. That is the reported
    # failure, and `equilibrate_certified` from that state now certifies at
    # 74.1899 cm3 with the magnesium back where it belongs. Too heavy for this
    # suite; what is asserted here is the mechanism, on a system of eight
    # species.
    sp4 = Dict(
        symbol(s) => s for s in build_species(
                datapath("slop98-inorganic-thermofun.json")
            )
    )
    cs4 = ChemicalSystem(
        [sp4[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
        ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
    )
    A4 = Float64.(cs4.SM.A)
    model4 = DiluteSolutionModel()

    # At a certified equilibrium every phase present is at LogSI = 0 and no
    # absent one is supersaturated. That is the built-in check on the whole
    # computation: if the present phases are not at zero, nothing else in the
    # result means anything.
    st4 = ChemicalState(cs4)
    set_quantity!(st4, "H2O@", 1.0u"kg")
    set_quantity!(st4, "Cal", 1.0e-3u"mol")
    eq4, cert4 = equilibrate_certified(st4; model = model4)
    @test cert4.optimal
    si4 = saturation_indices(eq4, model4)
    @test si4 isa AbstractDict
    @test length(si4) == length(cs4.species)
    @test abs(si4["Cal"]) < 1.0e-8            # present, hence saturated
    # Nothing to repair at an answer that certifies, and the route must say so
    # rather than invent a start.
    @test ChemistryLab._repair_start(
        eq4, model4, A4 * ustrip.(us"mol", eq4.n), 1.0e-16
    ) === nothing

    # Now a state where calcite is absent and the solution is grossly
    # supersaturated with respect to it — the shape of the reported failure,
    # without its size.
    st5 = ChemicalState(cs4)
    set_quantity!(st5, "H2O@", 1.0u"kg")
    set_quantity!(st5, "Ca+2", 0.1u"mol")
    set_quantity!(st5, "CO3-2", 0.1u"mol")
    si5 = saturation_indices(st5, model4)
    @test si5["Cal"] > 1.0
    b5 = A4 * ustrip.(us"mol", st5.n)

    fixed = ChemistryLab._repair_start(st5, model4, b5, 1.0e-16)
    @test fixed isa ChemicalState
    n5 = ustrip.(us"mol", fixed.n)
    i_cal = findfirst(s -> symbol(s) == "Cal", cs4.species)
    # The amount is what the recipe could make of it, scaled: a chemical bound,
    # not a guess at the answer. Calcite takes one Ca and one CO3, and there are
    # 0.1 mol of each.
    @test n5[i_cal] ≈ ChemistryLab._REPAIR_FRACTION * 0.1 rtol = 1.0e-8
    @test n5[i_cal] > 0
    # Everything else is left alone: the element balance of a start is not this
    # function's business, since `b` is fixed by the caller and every start is
    # projected onto it.
    @test all(
        n5[i] == ustrip(us"mol", st5.n[i]) for i in eachindex(n5) if i != i_cal
    )
    # And the temperature and pressure of the state it came from are carried.
    @test temperature(fixed) == temperature(st5)
    @test pressure(fixed) == pressure(st5)

end
