# Agreement with published values.
#
# These are not internal-consistency checks: each assertion pins a number that
# an outside source reports, so a drift in the database, the activity model or
# the solver shows up as a disagreement with the literature rather than as a
# silently changed figure in the documentation.

@testsection "Agreement with published values" begin

    @testsection "malonic acid pKa from SLOP98" begin
        # The titration example long called these species "maleic acid". They are
        # not: SLOP98 carries MALONIC-ACID,AQ / H-MALONATE,AQ / MALONATE,AQ, the
        # three-carbon diacid. The pKa derived from their ΔₐG⁰ settle the matter,
        # since malonic (2.83, 5.69) and maleic (1.92, 6.27) are far apart.
        sp = Dict(
            symbol(s) => s for s in vcat(
                    build_species(datapath("slop98-inorganic-thermofun.json")),
                    build_species(datapath("slop98-organic-thermofun.json"))
                )
        )

        for name in ("MalH2@", "MalH-", "Mal-2")
            @test haskey(sp, name)
        end

        # Three carbons, not four: malonate is C₃H₂O₄²⁻ at 102.05 g/mol, maleate
        # is C₄H₂O₄²⁻ at 114.06. The formula settles it on its own.
        @test ustrip(us"g/mol", sp["Mal-2"].M) ≈ 102.045 atol = 0.01
        @test ustrip(us"g/mol", sp["MalH2@"].M) ≈ 104.061 atol = 0.01
        @test occursin("H2O4", string(formula(sp["Mal-2"])))

        T = 298.15
        RT = 8.31446261815324 * T
        _pKa(num, den) = -log10(
            exp(
                -(sum(sp[x].ΔₐG⁰(T = T) for x in num) - sp[den].ΔₐG⁰(T = T)) / RT
            )
        )

        pKa1 = _pKa(("MalH-", "H+"), "MalH2@")
        pKa2 = _pKa(("Mal-2", "H+"), "MalH-")

        # Pinned to the database, tightly.
        @test pKa1 ≈ 2.851 atol = 0.002
        @test pKa2 ≈ 5.696 atol = 0.002

        # And agreeing with the tabulated values for malonic acid to 0.02,
        # while being nowhere near maleic acid's.
        @test abs(pKa1 - 2.83) < 0.03
        @test abs(pKa2 - 5.69) < 0.03
        @test abs(pKa1 - 1.92) > 0.5      # would hold if these were maleate
        @test abs(pKa2 - 6.27) > 0.5
    end

    @testsection "calcite: retrograde Kₛₚ in a closed system" begin
        # The documented 10-30 °C sweep once ran on `ChemicalState(cs)`, an empty
        # state: every element balance was zero, so all 21 points returned the
        # same trivial answer and the figure was a set of flat lines. This pins
        # the two things that must actually vary, and the one that must not.
        sp = Dict(
            symbol(s) => s for s in build_species(
                    datapath("slop98-inorganic-thermofun.json")
                )
        )
        species = [sp[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")]
        cs = ChemicalSystem(species, ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
        idx = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))

        function equilibrated(θ)
            st = ChemicalState(cs)
            set_quantity!(st, "Cal", 1.0e-3u"mol")
            set_quantity!(st, "H2O@", 1.0u"kg")
            V = volume(st)
            set_quantity!(st, "H+", 1.0e-4u"mol/L" * V.liquid)
            set_quantity!(st, "OH-", 1.0e-10u"mol/L" * V.liquid)
            set_temperature!(st, (273.15 + θ) * u"K")
            return equilibrate(st, OptimaOptimizer())
        end

        eq10, eq30 = equilibrated(10.0), equilibrated(30.0)

        # NOTE: one of these two solves reports `MaxIters` and still lands on the
        # right answer — the ionic product below matches the database Kₛₚ to 0.01
        # log units. The interior-point iteration stalls short of its own
        # stationarity test rather than short of the solution, which is the same
        # symptom the trace species showed in `equilibrium_reference.jl`. The
        # assertions here are on the composition, deliberately, so they measure
        # the answer and not the exit code.

        # The state is not empty: about a seventh of the calcite dissolves.
        @test ustrip(us"mol", eq10.n[idx["Ca+2"]]) > 1.0e-5

        molar(eq, s) = ustrip(us"mol", eq.n[idx[s]]) / ustrip(volume(eq).liquid) / 1000

        # Kₛₚ is retrograde: the ionic product the system settles on falls with T.
        logQ10 = log10(molar(eq10, "Ca+2") * molar(eq10, "CO3-2"))
        logQ30 = log10(molar(eq30, "Ca+2") * molar(eq30, "CO3-2"))
        @test logQ10 ≈ -8.411 atol = 0.01
        @test logQ30 ≈ -8.517 atol = 0.01
        @test logQ30 < logQ10 - 0.05

        # It equals the Kₛₚ the database itself gives for Cal = Ca+2 + CO3-2.
        for (θ, logQ) in ((10.0, logQ10), (30.0, logQ30))
            T = 273.15 + θ
            ΔG = (
                sp["Ca+2"].ΔₐG⁰(T = T) + sp["CO3-2"].ΔₐG⁰(T = T)
                    - sp["Cal"].ΔₐG⁰(T = T)
            )
            @test logQ ≈ -ΔG / (8.31446261815324 * T) / log(10) atol = 0.01
        end

        # And yet dissolved calcium RISES, because the pH falls and shifts
        # carbonate to bicarbonate. This is the counter-intuitive part the
        # documentation used to state backwards, so it is pinned by sign.
        @test pH(eq30) < pH(eq10) - 0.3
        @test molar(eq30, "Ca+2") > molar(eq10, "Ca+2")
        @test molar(eq30, "CO3-2") < molar(eq10, "CO3-2")
    end

end

@testsection "NaCl activity and osmotic coefficients against Hamer & Wu (1972)" begin
    # Hamer & Wu, *Osmotic Coefficients and Mean Activity Coefficients of
    # Uni-univalent Electrolytes in Water at 25 °C*, J. Phys. Chem. Ref. Data
    # 1(4), 1047-1100 (1972), Table 16 — a critical compilation, not a single
    # experiment. Read from the rendered page; the value at 6.144 mol/kg is
    # their saturated solution and is not used.
    #
    # This is the only assertion in the package that pins an activity model to a
    # *measurement* above a millimolal. Everything else about `PitzerActivityModel`
    # is internal consistency — the limiting law, Gibbs-Duhem — which cannot
    # distinguish a correct parameter set from a self-consistent wrong one. It
    # therefore also checks the transcription of Reardon's Na/Cl parameters:
    # nothing mistyped reproduces a measured curve over four decades.
    HAMER_WU_NACL = [
        # m [mol/kg]   φ        γ±
        (0.001, 0.988, 0.965),
        (0.01, 0.968, 0.903),
        (0.1, 0.933, 0.779),
        (0.5, 0.921, 0.681),
        (1.0, 0.936, 0.657),
        (2.0, 0.984, 0.668),
        (3.0, 1.045, 0.714),
        (4.0, 1.116, 0.783),
        (5.0, 1.191, 0.874),
        (6.0, 1.27, 0.986),
    ]

    subs = build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false)
    d = Dict(symbol(s) => s for s in subs)
    cs = ChemicalSystem([d[s] for s in split("H2O@ Na+ Cl-")], ["H2O@", "Na+", "Cl-"])
    M_w = ustrip(us"kg/mol", cs.species[only(cs.idx_solvent)][:M])
    n_w = 1 / M_w
    p = (ΔₐG⁰overRT = zeros(3), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)

    pz = activity_model(
        cs,
        PitzerActivityModel(;
            parameters = build_pitzer_parameters(datapath("pitzer-reardon1990.toml"))
        ),
    )
    bd = activity_model(cs, HKFActivityModel())

    worst_pz = 0.0
    for (m, φ_meas, γ_meas) in HAMER_WU_NACL
        out = pz([n_w, m, m], p)
        γ = exp((out[2] + out[3]) / 2 - log(m))
        φ = -out[1] / (M_w * 2m)
        # Both halves of the model, against measurement, independently.
        @test isapprox(γ, γ_meas; rtol = 0.01)
        @test isapprox(φ, φ_meas; rtol = 0.01)
        worst_pz = max(worst_pz, abs(γ - γ_meas) / γ_meas)
    end
    # Tighter than the assertions above: the whole curve, four decades of
    # molality, within half a percent — including the minimum near 1.2 mol/kg
    # and the return above unity at 6 mol/kg, neither of which a Debye-Hückel
    # form can produce at all.
    @test worst_pz < 0.005

    # The B-dot model on the same data. It is not being criticized for failing
    # outside its stated range; the point is that the range is real and that
    # nothing in its output announces the exit.
    dev(m, γ_meas) = let out = bd([n_w, m, m], p)
        abs(exp((out[2] + out[3]) / 2 - log(m)) - γ_meas) / γ_meas
    end
    @test dev(0.001, 0.965) < 0.01        # agrees where it should
    @test dev(0.1, 0.779) > 0.03        # 5 % out at a tenth molal
    @test dev(1.0, 0.657) > 0.15        # 19 % at one
    @test dev(6.0, 0.986) > 0.35        # 44 % at six
end
