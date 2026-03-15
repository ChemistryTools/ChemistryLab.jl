@testsection "SolidSolution" begin

    # ── Shared fixtures ───────────────────────────────────────────────────────
    h2o       = Species("H2O";       aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    calcite   = Species("CaCO3";     aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    magnesite = Species("MgCO3";     aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    dolomite  = Species("CaMg(CO3)2"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)

    R_gas = ustrip(us"J/mol/K", Constants.R)  # ≈ 8.31446 J/(mol·K)
    T25   = 298.15                          # K
    RT25  = R_gas * T25

    # ── IdealSolidSolutionMixingModel ─────────────────────────────────────────
    @testsection "IdealSolidSolutionMixingModel" begin
        m = IdealSolidSolutionMixingModel()
        @test m isa IdealSolidSolutionMixingModel
        @test lnγ_ss(m, [0.7, 0.3], RT25) == [0.0, 0.0]
        @test lnγ_ss(m, [1.0 / 3, 1.0 / 3, 1.0 / 3], RT25) == [0.0, 0.0, 0.0]
    end

    # ── RegularSolidSolutionMixingModel — constructor ─────────────────────────
    @testsection "RegularSolidSolutionMixingModel construction" begin
        W_ok = [0.0 3000.0; 3000.0 0.0]
        @test RegularSolidSolutionMixingModel(W_ok) isa RegularSolidSolutionMixingModel

        # Asymmetric W
        @test_throws ArgumentError RegularSolidSolutionMixingModel(
            [0.0 1000.0; 2000.0 0.0],
        )
        # Non-square W
        @test_throws ArgumentError RegularSolidSolutionMixingModel(
            [0.0 1000.0 500.0; 1000.0 0.0 800.0],
        )
        # Non-zero diagonal
        @test_throws ArgumentError RegularSolidSolutionMixingModel(
            [1.0 3000.0; 3000.0 0.0],
        )
    end

    # ── RegularSolidSolutionMixingModel — lnγ_ss values ──────────────────────
    @testsection "RegularSolidSolutionMixingModel lnγ_ss" begin
        W   = 3000.0
        mix = RegularSolidSolutionMixingModel([0.0 W; W 0.0])

        # Binary equimolar: RT ln γᵢ = W × 0.5² = W / 4
        x = [0.5, 0.5]
        lnγ = lnγ_ss(mix, x, RT25)
        @test lnγ[1] ≈ W * 0.5^2 / RT25  rtol = 1e-10
        @test lnγ[2] ≈ W * 0.5^2 / RT25  rtol = 1e-10

        # Binary asymmetric: RT ln γ₁ = W x₂²
        x2 = [0.3, 0.7]
        lnγ2 = lnγ_ss(mix, x2, RT25)
        @test lnγ2[1] ≈ W * 0.7^2 / RT25  rtol = 1e-10
        @test lnγ2[2] ≈ W * 0.3^2 / RT25  rtol = 1e-10

        # Dilute limit: x₁ → 1 ⟹ ln γ₁ → 0, ln γ₂ → W/RT (Henry's law)
        x_dilute = [1.0 - 1e-10, 1e-10]
        lnγ_dilute = lnγ_ss(mix, x_dilute, RT25)
        @test abs(lnγ_dilute[1]) < 1e-6
        @test lnγ_dilute[2] ≈ W / RT25  rtol = 1e-6

        # Gibbs–Duhem consistency: x₁ Δ(ln γ₁)/Δx₁ + x₂ Δ(ln γ₂)/Δx₁ = 0
        # For symmetric Margules: d(ln γ₁)/dx₁ = -2W x₂/RT, d(ln γ₂)/dx₁ = 2W x₁/RT
        # ⟹ x₁(-2W x₂) + x₂(2W x₁) = 0 ✓  (verified analytically)
    end

    # ── SolidSolution constructor ─────────────────────────────────────────────
    @testsection "SolidSolution constructor" begin
        # Valid ideal binary
        ss = SolidSolution(["CaCO3", "MgCO3"], IdealSolidSolutionMixingModel())
        @test ss isa SolidSolution
        @test length(ss.end_member_symbols) == 2

        # Empty symbol list
        @test_throws ArgumentError SolidSolution(
            String[], IdealSolidSolutionMixingModel(),
        )

        # W size mismatch (3×3 for 2 end-members)
        W_bad = zeros(3, 3)
        @test_throws DimensionMismatch SolidSolution(
            ["CaCO3", "MgCO3"], RegularSolidSolutionMixingModel(W_bad),
        )
    end

    # ── SolidSolutionActivityModel — validation at activity_model time ────────
    @testsection "SolidSolutionActivityModel validation" begin
        cs = ChemicalSystem([h2o, calcite, magnesite])

        # Valid: all symbols present and AS_CRYSTAL
        ss    = SolidSolution(["CaCO3", "MgCO3"], IdealSolidSolutionMixingModel())
        model = SolidSolutionActivityModel(DiluteSolutionModel(), [ss])
        @test activity_model(cs, model) isa Function

        # Symbol not in ChemicalSystem
        ss_bad = SolidSolution(["CaCO3", "NonExistent"], IdealSolidSolutionMixingModel())
        @test_throws ArgumentError activity_model(
            cs, SolidSolutionActivityModel(DiluteSolutionModel(), [ss_bad]),
        )

        # Aqueous species used as end-member
        ss_aq = SolidSolution(["CaCO3", "H2O"], IdealSolidSolutionMixingModel())
        @test_throws ArgumentError activity_model(
            cs, SolidSolutionActivityModel(DiluteSolutionModel(), [ss_aq]),
        )

        # Same end-member in two groups
        ss1 = SolidSolution(["CaCO3", "MgCO3"], IdealSolidSolutionMixingModel())
        ss2 = SolidSolution(["CaCO3"], IdealSolidSolutionMixingModel())
        @test_throws ArgumentError activity_model(
            cs, SolidSolutionActivityModel(DiluteSolutionModel(), [ss1, ss2]),
        )
    end

    # ── Ideal binary SS — activity values ─────────────────────────────────────
    @testsection "Ideal binary SS activities" begin
        cs  = ChemicalSystem([h2o, calcite, magnesite])
        ss  = SolidSolution(["CaCO3", "MgCO3"], IdealSolidSolutionMixingModel())
        lna = activity_model(cs, SolidSolutionActivityModel(DiluteSolutionModel(), [ss]))

        # x_calcite = 0.75, x_magnesite = 0.25
        n = [55.5, 0.75, 0.25]
        p = (ΔₐG⁰overT = zeros(3), ϵ = 1e-30, T = T25)
        out = lna(n, p)

        i_cal = findfirst(s -> symbol(s) == "CaCO3",  cs.species)
        i_mag = findfirst(s -> symbol(s) == "MgCO3",  cs.species)

        @test out[i_cal] ≈ log(0.75)  rtol = 1e-10
        @test out[i_mag] ≈ log(0.25)  rtol = 1e-10
    end

    # ── Pure crystal not in any SS keeps ln a = 0 ─────────────────────────────
    @testsection "Pure crystal outside SS unaffected" begin
        cs  = ChemicalSystem([h2o, calcite, magnesite, dolomite])
        ss  = SolidSolution(["CaCO3", "MgCO3"], IdealSolidSolutionMixingModel())
        lna = activity_model(cs, SolidSolutionActivityModel(DiluteSolutionModel(), [ss]))

        n = [55.5, 0.5, 0.5, 0.1]
        p = (ΔₐG⁰overT = zeros(4), ϵ = 1e-30, T = T25)
        out = lna(n, p)

        i_dol = findfirst(s -> symbol(s) == "CaMg(CO3)2", cs.species)
        @test out[i_dol] == 0.0
    end

    # ── SS entirely absent: no crash, ln a = 0 ────────────────────────────────
    @testsection "SS entirely absent" begin
        cs  = ChemicalSystem([h2o, calcite, magnesite])
        ss  = SolidSolution(["CaCO3", "MgCO3"], IdealSolidSolutionMixingModel())
        lna = activity_model(cs, SolidSolutionActivityModel(DiluteSolutionModel(), [ss]))

        # Both crystals absent (well below ϵ threshold)
        n = [55.5, 1e-35, 1e-35]
        p = (ΔₐG⁰overT = zeros(3), ϵ = 1e-30, T = T25)
        out = lna(n, p)

        i_cal = findfirst(s -> symbol(s) == "CaCO3", cs.species)
        i_mag = findfirst(s -> symbol(s) == "MgCO3", cs.species)
        @test out[i_cal] == 0.0
        @test out[i_mag] == 0.0
    end

    # ── Single end-member SS (x = 1) → ln a = 0 ──────────────────────────────
    @testsection "Single end-member SS" begin
        cs  = ChemicalSystem([h2o, calcite])
        ss  = SolidSolution(["CaCO3"], IdealSolidSolutionMixingModel())
        lna = activity_model(cs, SolidSolutionActivityModel(DiluteSolutionModel(), [ss]))

        n = [55.5, 1.0]
        p = (ΔₐG⁰overT = zeros(2), ϵ = 1e-30, T = T25)
        out = lna(n, p)

        i_cal = findfirst(s -> symbol(s) == "CaCO3", cs.species)
        @test out[i_cal] ≈ 0.0  atol = 1e-14
    end

    # ── Ideal ternary SS ──────────────────────────────────────────────────────
    @testsection "Ideal ternary SS" begin
        cs  = ChemicalSystem([h2o, calcite, magnesite, dolomite])
        ss  = SolidSolution(
            ["CaCO3", "MgCO3", "CaMg(CO3)2"], IdealSolidSolutionMixingModel(),
        )
        lna = activity_model(cs, SolidSolutionActivityModel(DiluteSolutionModel(), [ss]))

        # n_SS = 1 + 2 + 3 = 6 → x = [1/6, 2/6, 3/6]
        n = [55.5, 1.0, 2.0, 3.0]
        p = (ΔₐG⁰overT = zeros(4), ϵ = 1e-30, T = T25)
        out = lna(n, p)

        i_cal = findfirst(s -> symbol(s) == "CaCO3",      cs.species)
        i_mag = findfirst(s -> symbol(s) == "MgCO3",      cs.species)
        i_dol = findfirst(s -> symbol(s) == "CaMg(CO3)2", cs.species)

        @test out[i_cal] ≈ log(1.0 / 6)  rtol = 1e-10
        @test out[i_mag] ≈ log(2.0 / 6)  rtol = 1e-10
        @test out[i_dol] ≈ log(3.0 / 6)  rtol = 1e-10
    end

    # ── Regular binary SS — activity values ───────────────────────────────────
    @testsection "Regular binary SS activities" begin
        W     = 3000.0
        mix   = RegularSolidSolutionMixingModel([0.0 W; W 0.0])
        cs    = ChemicalSystem([h2o, calcite, magnesite])
        ss    = SolidSolution(["CaCO3", "MgCO3"], mix)
        lna   = activity_model(cs, SolidSolutionActivityModel(DiluteSolutionModel(), [ss]))

        # x_calcite = 0.3, x_magnesite = 0.7
        n = [55.5, 0.3, 0.7]
        p = (ΔₐG⁰overT = zeros(3), ϵ = 1e-30, T = T25)
        out = lna(n, p)

        i_cal = findfirst(s -> symbol(s) == "CaCO3", cs.species)
        i_mag = findfirst(s -> symbol(s) == "MgCO3", cs.species)

        # ln aᵢ = ln xᵢ + W xⱼ² / RT
        RT = R_gas * T25
        @test out[i_cal] ≈ log(0.3) + W * 0.7^2 / RT  rtol = 1e-10
        @test out[i_mag] ≈ log(0.7) + W * 0.3^2 / RT  rtol = 1e-10
    end

    # ── Two independent SS groups ─────────────────────────────────────────────
    @testsection "Two independent SS groups" begin
        iron_ox = Species("Fe2O3"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
        cs      = ChemicalSystem([h2o, calcite, magnesite, iron_ox])

        ss1 = SolidSolution(["CaCO3", "MgCO3"], IdealSolidSolutionMixingModel())
        ss2 = SolidSolution(["Fe2O3"], IdealSolidSolutionMixingModel())
        lna = activity_model(
            cs, SolidSolutionActivityModel(DiluteSolutionModel(), [ss1, ss2]),
        )

        n = [55.5, 1.0, 1.0, 2.0]
        p = (ΔₐG⁰overT = zeros(4), ϵ = 1e-30, T = T25)
        out = lna(n, p)

        i_cal  = findfirst(s -> symbol(s) == "CaCO3",  cs.species)
        i_mag  = findfirst(s -> symbol(s) == "MgCO3",  cs.species)
        i_iron = findfirst(s -> symbol(s) == "Fe2O3",  cs.species)

        @test out[i_cal]  ≈ log(0.5)  rtol = 1e-10
        @test out[i_mag]  ≈ log(0.5)  rtol = 1e-10
        @test out[i_iron] ≈ 0.0  atol = 1e-14  # singleton SS, x = 1
    end

    # ── Regular ternary lnγ_ss — equimolar and asymmetric ─────────────────────
    @testsection "Regular ternary SS lnγ_ss" begin
        W     = 6000.0
        W_mat = [0.0 W W; W 0.0 W; W W 0.0]
        mix   = RegularSolidSolutionMixingModel(W_mat)

        # Equimolar x = [1/3, 1/3, 1/3]:
        # RT ln γᵢ = W(1/3)(2/3) + W(1/3)(2/3) - W(1/3)(1/3) = W/3
        x_eq = [1.0 / 3, 1.0 / 3, 1.0 / 3]
        lnγ_eq = lnγ_ss(mix, x_eq, RT25)
        @test lnγ_eq[1] ≈ W / 3 / RT25  rtol = 1e-10
        @test lnγ_eq[2] ≈ W / 3 / RT25  rtol = 1e-10
        @test lnγ_eq[3] ≈ W / 3 / RT25  rtol = 1e-10

        # Asymmetric x = [0.5, 0.3, 0.2] (same W everywhere):
        # RT ln γ₁ = W x₂(1-x₁) + W x₃(1-x₁) - W x₂x₃
        #          = W(0.3·0.5 + 0.2·0.5 − 0.3·0.2) = 0.19 W
        # RT ln γ₂ = W x₁(1-x₂) + W x₃(1-x₂) - W x₁x₃
        #          = W(0.5·0.7 + 0.2·0.7 − 0.5·0.2) = 0.39 W
        # RT ln γ₃ = W x₁(1-x₃) + W x₂(1-x₃) - W x₁x₂
        #          = W(0.5·0.8 + 0.3·0.8 − 0.5·0.3) = 0.49 W
        x_as = [0.5, 0.3, 0.2]
        lnγ_as = lnγ_ss(mix, x_as, RT25)
        @test lnγ_as[1] ≈ W * (0.3 * 0.5 + 0.2 * 0.5 - 0.3 * 0.2) / RT25  rtol = 1e-10
        @test lnγ_as[2] ≈ W * (0.5 * 0.7 + 0.2 * 0.7 - 0.5 * 0.2) / RT25  rtol = 1e-10
        @test lnγ_as[3] ≈ W * (0.5 * 0.8 + 0.3 * 0.8 - 0.5 * 0.3) / RT25  rtol = 1e-10
    end

    # ── Gibbs–Duhem numerical check (binary regular) ──────────────────────────
    # At constant T,P: Σᵢ xᵢ d(ln γᵢ)/dx₁ = 0.
    # For symmetric Margules: x₁(-2Wx₂/RT) + x₂(+2Wx₁/RT) = 0 exactly.
    @testsection "Gibbs-Duhem consistency (binary regular)" begin
        W   = 5000.0
        mix = RegularSolidSolutionMixingModel([0.0 W; W 0.0])
        Δ   = 1e-7

        for x1 in [0.1, 0.3, 0.5, 0.7, 0.9]
            x2 = 1 - x1
            lnγ_fwd = lnγ_ss(mix, [x1 + Δ, x2 - Δ], RT25)
            lnγ_bwd = lnγ_ss(mix, [x1 - Δ, x2 + Δ], RT25)
            d_lnγ1 = (lnγ_fwd[1] - lnγ_bwd[1]) / (2Δ)
            d_lnγ2 = (lnγ_fwd[2] - lnγ_bwd[2]) / (2Δ)
            @test abs(x1 * d_lnγ1 + x2 * d_lnγ2) < 1e-5
        end
    end

    # ── W = 0 regular solution equals ideal ───────────────────────────────────
    @testsection "Regular with W=0 equals ideal" begin
        mix_zero  = RegularSolidSolutionMixingModel([0.0 0.0; 0.0 0.0])
        mix_ideal = IdealSolidSolutionMixingModel()

        for x in [[0.3, 0.7], [0.5, 0.5], [0.9, 0.1]]
            @test lnγ_ss(mix_zero, x, RT25) == lnγ_ss(mix_ideal, x, RT25)
        end
    end

    # ── Regular ternary SS — full activity_model closure ──────────────────────
    @testsection "Regular ternary SS in activity_model" begin
        W     = 3000.0
        W_mat = [0.0 W W; W 0.0 W; W W 0.0]
        mix   = RegularSolidSolutionMixingModel(W_mat)
        cs    = ChemicalSystem([h2o, calcite, magnesite, dolomite])
        ss    = SolidSolution(["CaCO3", "MgCO3", "CaMg(CO3)2"], mix)
        lna   = activity_model(cs, SolidSolutionActivityModel(DiluteSolutionModel(), [ss]))

        # n_SS = 5 + 3 + 2 = 10 → x = [0.5, 0.3, 0.2]
        n = [55.5, 5.0, 3.0, 2.0]
        p = (ΔₐG⁰overT = zeros(4), ϵ = 1e-30, T = T25)
        out = lna(n, p)

        RT = R_gas * T25
        i_cal = findfirst(s -> symbol(s) == "CaCO3",      cs.species)
        i_mag = findfirst(s -> symbol(s) == "MgCO3",      cs.species)
        i_dol = findfirst(s -> symbol(s) == "CaMg(CO3)2", cs.species)

        x1, x2, x3 = 0.5, 0.3, 0.2
        lnγ1 = W * (x2 * (1 - x1) + x3 * (1 - x1) - x2 * x3) / RT
        lnγ2 = W * (x1 * (1 - x2) + x3 * (1 - x2) - x1 * x3) / RT
        lnγ3 = W * (x1 * (1 - x3) + x2 * (1 - x3) - x1 * x2) / RT

        @test out[i_cal] ≈ log(x1) + lnγ1  rtol = 1e-10
        @test out[i_mag] ≈ log(x2) + lnγ2  rtol = 1e-10
        @test out[i_dol] ≈ log(x3) + lnγ3  rtol = 1e-10
    end

end  # @testsection "SolidSolution"
