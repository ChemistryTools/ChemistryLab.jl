using ChemistryLab
using DynamicQuantities
using ForwardDiff
using LinearAlgebra
using Test

# ── Helpers ────────────────────────────────────────────────────────────────────

# Build a minimal Species for solid solution end-member tests (no thermodynamic data needed)
_ss_em(sym) = Species(sym; aggregate_state = AS_CRYSTAL, class = SC_SSENDMEMBER)

# ── Tests ──────────────────────────────────────────────────────────────────────

@testsection "SC_SSENDMEMBER enum" begin
    @test SC_SSENDMEMBER isa Class
    @test SC_SSENDMEMBER != SC_COMPONENT
    @test SC_SSENDMEMBER != SC_AQSOLUTE
    @test SC_SSENDMEMBER != SC_UNDEF
    sp = _ss_em("AFm1")
    @test class(sp) == SC_SSENDMEMBER
    @test aggregate_state(sp) == AS_CRYSTAL
end

@testsection "IdealSolidSolutionModel constructor" begin
    m = IdealSolidSolutionModel()
    @test m isa IdealSolidSolutionModel
    @test m isa AbstractSolidSolutionModel
end

@testsection "RedlichKisterModel constructor" begin
    # Default all-zero
    m0 = RedlichKisterModel()
    @test m0 isa RedlichKisterModel{Float64}
    @test m0.a0 == 0.0 && m0.a1 == 0.0 && m0.a2 == 0.0

    # Custom parameters
    m1 = RedlichKisterModel(a0 = 4000.0, a1 = 500.0)
    @test m1.a0 ≈ 4000.0
    @test m1.a1 ≈ 500.0
    @test m1.a2 == 0.0

    # Float32 promotion (all three args must be Float32)
    m32 = RedlichKisterModel(a0 = 1000.0f0, a1 = 200.0f0, a2 = 50.0f0)
    @test m32 isa RedlichKisterModel{Float32}

    # Mixed promotes to Float64
    m_mixed = RedlichKisterModel(a0 = 1000.0f0, a1 = 200.0)
    @test m_mixed isa RedlichKisterModel{Float64}
end

@testsection "SolidSolutionPhase constructor" begin
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")

    # Ideal (default)
    ss = SolidSolutionPhase("SS", [em1, em2])
    @test ss isa SolidSolutionPhase
    @test name(ss) == "SS"
    @test length(end_members(ss)) == 2
    @test model(ss) isa IdealSolidSolutionModel

    # Explicit model
    rk = RedlichKisterModel(a0 = 3000.0)
    ss_rk = SolidSolutionPhase("SS_RK", [em1, em2]; model = rk)
    @test model(ss_rk) isa RedlichKisterModel

    # Error: wrong aggregate_state
    em_bad = Species("Bad"; aggregate_state = AS_AQUEOUS, class = SC_SSENDMEMBER)
    @test_throws ErrorException SolidSolutionPhase("Bad", [em_bad])

    # Auto-requalification: SC_COMPONENT end-members are silently promoted to SC_SSENDMEMBER
    em_comp = Species("Comp"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    ss_auto = @test_nowarn SolidSolutionPhase("Auto", [em_comp])
    @test class(end_members(ss_auto)[1]) == SC_SSENDMEMBER
    @test class(em_comp) == SC_COMPONENT   # original unchanged

    # Error: Redlich-Kister with != 2 end-members
    em3 = _ss_em("Em3")
    @test_throws ErrorException SolidSolutionPhase(
        "SS3", [em1, em2, em3]; model = RedlichKisterModel()
    )
    @test_throws ErrorException SolidSolutionPhase(
        "SS1", [em1]; model = RedlichKisterModel()
    )
end

@testsection "ChemicalSystem without solid solutions" begin
    em1 = _ss_em("Em1")
    cs = ChemicalSystem([em1])
    # Backward compatibility — solid_solutions field is Nothing
    @test cs.solid_solutions === nothing
    @test isempty(cs.ss_groups)
    @test isempty(cs.idx_ssendmembers)
end

@testsection "ChemicalSystem with solid solutions" begin
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    em3 = _ss_em("Em3")
    ss1 = SolidSolutionPhase("SS1", [em1, em2])
    ss2 = SolidSolutionPhase("SS2", [em3])

    cs = ChemicalSystem([em1, em2, em3]; solid_solutions = [ss1, ss2])

    @test !isnothing(cs.solid_solutions)
    @test length(cs.solid_solutions) == 2
    @test cs.ss_groups == [[1, 2], [3]]
    @test cs.idx_ssendmembers == [1, 2, 3]

    # All end-members are in idx_crystal (AS_CRYSTAL aggregate state)
    @test Set(cs.idx_ssendmembers) ⊆ Set(cs.idx_crystal)

    # Accessor
    @test solid_solutions(cs) === cs.solid_solutions

    # Error: end-member not in species list
    em_extra = _ss_em("Extra")
    ss_bad = SolidSolutionPhase("Bad", [em_extra])
    @test_throws ErrorException ChemicalSystem([em1, em2]; solid_solutions = [ss_bad])
end

@testsection "Ideal SS binary: ln aᵢ = ln xᵢ" begin
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    ss = SolidSolutionPhase("SS", [em1, em2])
    cs = ChemicalSystem([em1, em2]; solid_solutions = [ss])

    lna = activity_model(cs, DiluteSolutionModel())
    p = (ΔₐG⁰overT = zeros(2), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)

    # x₁ = 0.7, x₂ = 0.3
    n = [0.7, 0.3]
    out = lna(n, p)
    @test out[1] ≈ log(0.7) rtol = 1.0e-8
    @test out[2] ≈ log(0.3) rtol = 1.0e-8

    # x₁ = 0.5, x₂ = 0.5 (symmetric)
    n2 = [0.5, 0.5]
    out2 = lna(n2, p)
    @test out2[1] ≈ log(0.5) rtol = 1.0e-8
    @test out2[2] ≈ log(0.5) rtol = 1.0e-8
    @test out2[1] ≈ out2[2] rtol = 1.0e-10
end

@testsection "Ideal SS: pure end-member (x₁ → 1)" begin
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    ss = SolidSolutionPhase("SS", [em1, em2])
    cs = ChemicalSystem([em1, em2]; solid_solutions = [ss])

    lna = activity_model(cs, DiluteSolutionModel())
    p = (ΔₐG⁰overT = zeros(2), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)

    # n₂ very small → x₁ ≈ 1, ln a₁ ≈ 0
    n_pure = [1.0 - 1.0e-12, 1.0e-12]
    out = lna(n_pure, p)
    @test out[1] ≈ 0.0 atol = 1.0e-8

    # n₁ very small → x₂ ≈ 1, ln a₂ ≈ 0
    n_pure2 = [1.0e-12, 1.0 - 1.0e-12]
    out2 = lna(n_pure2, p)
    @test out2[2] ≈ 0.0 atol = 1.0e-8
end

@testsection "Redlich-Kister binary: analytical values" begin
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    a0 = 4000.0  # J/mol
    a1 = 800.0   # J/mol
    rk = RedlichKisterModel(a0 = a0, a1 = a1)
    ss = SolidSolutionPhase("SS_RK", [em1, em2]; model = rk)
    cs = ChemicalSystem([em1, em2]; solid_solutions = [ss])

    T = 298.15
    RT = 8.31446261815324 * T

    lna = activity_model(cs, DiluteSolutionModel())
    p = (ΔₐG⁰overT = zeros(2), T = T, P = 1.0e5, ϵ = 1.0e-30)

    # x₁ = 0.3, x₂ = 0.7
    x1, x2 = 0.3, 0.7
    n = [x1, x2]
    out = lna(n, p)

    # Analytical Redlich-Kister: ln γ₁
    expected_lng1 = x2^2 * (a0 + a1 * (3 * x1 - x2)) / RT
    expected_lng2 = x1^2 * (a0 - a1 * (3 * x2 - x1)) / RT
    @test out[1] ≈ log(x1) + expected_lng1 rtol = 1.0e-8
    @test out[2] ≈ log(x2) + expected_lng2 rtol = 1.0e-8

    # At x = 0.5 (symmetric), a1=0 case
    rk0 = RedlichKisterModel(a0 = a0)
    ss0 = SolidSolutionPhase("SS0", [em1, em2]; model = rk0)
    cs0 = ChemicalSystem([em1, em2]; solid_solutions = [ss0])
    lna0 = activity_model(cs0, DiluteSolutionModel())
    out0 = lna0([0.5, 0.5], p)
    @test out0[1] ≈ out0[2] rtol = 1.0e-10   # symmetric at x=0.5 with a1=0
end

@testsection "SS in DiluteSolutionModel (mixed aqueous + SS)" begin
    H2O = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    NaCl = Species("NaCl"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    ss = SolidSolutionPhase("SS", [em1, em2])

    cs = ChemicalSystem([H2O, NaCl, em1, em2]; solid_solutions = [ss])
    lna = activity_model(cs, DiluteSolutionModel())

    n_w = 55.5
    n_Na = 0.1
    n1, n2 = 0.6, 0.4
    n = [n_w, n_Na, n1, n2]
    p = (ΔₐG⁰overT = zeros(4), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    out = lna(n, p)

    # SS end-members: should match ideal mixing formula
    @test out[3] ≈ log(n1 / (n1 + n2)) rtol = 1.0e-8
    @test out[4] ≈ log(n2 / (n1 + n2)) rtol = 1.0e-8

    # Aqueous species: should be non-zero and unaffected by SS
    @test out[1] != 0.0   # solvent
    @test out[2] != 0.0   # solute
end

@testsection "SS in HKFActivityModel (mixed aqueous + SS)" begin
    H2O = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    ss = SolidSolutionPhase("SS", [em1, em2])

    cs = ChemicalSystem([H2O, em1, em2]; solid_solutions = [ss])
    lna = activity_model(cs, HKFActivityModel())

    n_w = 55.5
    n1, n2 = 0.6, 0.4
    n = [n_w, n1, n2]
    p = (ΔₐG⁰overT = zeros(3), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    out = lna(n, p)

    # SS end-members correct
    @test out[2] ≈ log(n1 / (n1 + n2)) rtol = 1.0e-8
    @test out[3] ≈ log(n2 / (n1 + n2)) rtol = 1.0e-8
end

@testsection "Gibbs-Duhem local SS (ideal)" begin
    # For an ideal SS phase, Σᵢ xᵢ dln(aᵢ)/dξ = 0 where ξ shifts composition
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    ss = SolidSolutionPhase("SS", [em1, em2])
    cs = ChemicalSystem([em1, em2]; solid_solutions = [ss])

    lna = activity_model(cs, DiluteSolutionModel())
    p = (ΔₐG⁰overT = zeros(2), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)

    n0 = [0.7, 0.3]
    # Perturbation: shift n₁ ↑ ε, n₂ ↓ ε (internal composition change, total constant)
    ξ = 1.0e-6
    n⁺ = [0.7 + ξ, 0.3 - ξ]
    n⁻ = [0.7 - ξ, 0.3 + ξ]
    dμ = (lna(n⁺, p) - lna(n⁻, p)) / (2 * ξ)   # finite difference

    # Gibbs-Duhem: Σ nᵢ dμᵢ/dξ ≈ 0 (exact for ideal SS)
    residual = abs(sum(n0 .* dμ)) / max(norm(n0 .* abs.(dμ)), 1.0)
    @test residual < 1.0e-6
end

@testsection "ForwardDiff gradient (SS + aqueous)" begin
    H2O = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    Na = Species("Na+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    ss = SolidSolutionPhase("SS", [em1, em2])
    cs = ChemicalSystem([H2O, Na, em1, em2]; solid_solutions = [ss])

    lna = activity_model(cs, DiluteSolutionModel())
    p = (ΔₐG⁰overT = zeros(4), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    n0 = [55.5, 0.1, 0.7, 0.3]

    # Gradient of each component
    for i in 1:4
        g = ForwardDiff.gradient(n -> lna(n, p)[i], n0)
        @test all(isfinite, g)
    end

    # Jacobian of the full lna vector
    J = ForwardDiff.jacobian(n -> lna(n, p), n0)
    @test all(isfinite, J)
end

@testsection "ForwardDiff gradient (pure SS, no aqueous)" begin
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    ss = SolidSolutionPhase("SS", [em1, em2])
    cs = ChemicalSystem([em1, em2]; solid_solutions = [ss])

    lna = activity_model(cs, DiluteSolutionModel())
    p = (ΔₐG⁰overT = zeros(2), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    n0 = [0.7, 0.3]

    J = ForwardDiff.jacobian(n -> lna(n, p), n0)
    @test all(isfinite, J)
end

@testsection "ForwardDiff gradient (Redlich-Kister)" begin
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    rk = RedlichKisterModel(a0 = 3000.0, a1 = 500.0, a2 = 100.0)
    ss = SolidSolutionPhase("SS", [em1, em2]; model = rk)
    cs = ChemicalSystem([em1, em2]; solid_solutions = [ss])

    lna = activity_model(cs, DiluteSolutionModel())
    p = (ΔₐG⁰overT = zeros(2), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    n0 = [0.6, 0.4]

    J = ForwardDiff.jacobian(n -> lna(n, p), n0)
    @test all(isfinite, J)
end

@testsection "build_potentials with SS" begin
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    ss = SolidSolutionPhase("SS", [em1, em2])
    cs = ChemicalSystem([em1, em2]; solid_solutions = [ss])

    μ = build_potentials(cs, DiluteSolutionModel())
    n = [0.7, 0.3]
    ΔaGoT = [-100.0, -110.0]
    p = (ΔₐG⁰overT = ΔaGoT, T = 298.15, P = 1.0e5, ϵ = 1.0e-30)

    out = μ(n, p)
    # μᵢ/RT = ΔₐG⁰ᵢ/RT + ln aᵢ
    lna_out = activity_model(cs, DiluteSolutionModel())(n, p)
    @test out ≈ ΔaGoT .+ lna_out rtol = 1.0e-10
end

@testsection "RegularSolutionModel" begin
    # Construction and its guards.
    m = RegularSolutionModel([0.0 4000.0; 4000.0 0.0])
    @test m isa RegularSolutionModel
    @test m.W[1, 2] == 4000.0
    @test_throws ArgumentError RegularSolutionModel([0.0 1.0 2.0; 1.0 0.0 3.0])
    @test_throws ArgumentError RegularSolutionModel([0.0 4000.0; 1.0 0.0])

    # The binary limit is exactly Redlich-Kister with a0 = W and a1 = a2 = 0.
    # That is what makes this the right generalization rather than a new model:
    # `ln γ₁ = W x₂²/RT`, `ln γ₂ = W x₁²/RT`.
    T = 298.15
    RT = 8.31446261815324 * T
    W = 4000.0
    reg = RegularSolutionModel([0.0 W; W 0.0])
    rk = RedlichKisterModel(a0 = W)
    for x1 in (0.05, 0.3, 0.5, 0.7, 0.95)
        x = [x1, 1 - x1]
        for k in 1:2
            @test ChemistryLab._excess_ln_gamma(reg, k, x, T) ≈
                ChemistryLab._excess_ln_gamma(rk, k, x, T) rtol = 1.0e-12
        end
        @test ChemistryLab._excess_ln_gamma(reg, 1, x, T) ≈ W * x[2]^2 / RT rtol = 1.0e-12
        @test ChemistryLab._excess_ln_gamma(reg, 2, x, T) ≈ W * x[1]^2 / RT rtol = 1.0e-12
    end

    # W = 0 is ideal mixing, at any arity.
    ideal3 = RegularSolutionModel(zeros(3, 3))
    for k in 1:3
        @test ChemistryLab._excess_ln_gamma(ideal3, k, [0.2, 0.3, 0.5], T) ≈ 0.0 atol = 1.0e-15
    end

    # Ternary: the closed form, checked against the definition
    # `ln γ_k = Σ_{j≠k} W_kj x_j − Σ_{i<j} W_ij x_i x_j`, all over RT.
    W3 = [0.0 3000.0 -1500.0; 3000.0 0.0 2000.0; -1500.0 2000.0 0.0]
    reg3 = RegularSolutionModel(W3)
    x = [0.2, 0.3, 0.5]
    quad = sum(W3[i, j] * x[i] * x[j] for i in 1:3 for j in (i + 1):3) / RT
    for k in 1:3
        lin = sum(W3[k, j] * x[j] for j in 1:3 if j != k) / RT
        @test ChemistryLab._excess_ln_gamma(reg3, k, x, T) ≈ lin - quad rtol = 1.0e-12
    end

    # Gibbs-Duhem at fixed T: Σ x_k d(ln γ_k) = 0, so Σ x_k ln γ_k must equal
    # G^ex/RT. This is the property that makes the partial derivatives mutually
    # consistent, and it is what a hand-written `ln γ` most often gets wrong.
    gex = sum(W3[i, j] * x[i] * x[j] for i in 1:3 for j in (i + 1):3) / RT
    @test sum(x[k] * ChemistryLab._excess_ln_gamma(reg3, k, x, T) for k in 1:3) ≈
        gex rtol = 1.0e-12

    # A phase of any arity accepts it, unlike Redlich-Kister which is binary.
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    dict = Dict(symbol(s) => s for s in substances)
    six = ["CSHQ-TobD", "CSHQ-TobH", "CSHQ-JenH", "CSHQ-JenD", "KSiOH", "NaSiOH"]
    ss = SolidSolutionPhase(
        "CSHQ", [dict[m] for m in six]; model = RegularSolutionModel(zeros(6, 6))
    )
    @test length(end_members(ss)) == 6
    @test model(ss) isa RegularSolutionModel
    @test_throws ErrorException SolidSolutionPhase(
        "CSHQ", [dict[m] for m in six]; model = RedlichKisterModel(a0 = 1.0)
    )

    # And the shipped file loads, with the five phases added in 0.15.0.
    ss_all = build_solid_solutions(datapath("solid_solutions.toml"), dict)
    names = Set(p.name for p in ss_all)
    @test length(ss_all) == 11
    for n in (
            "Straetlingite_ss", "AFm_SO4_OH", "AFt_SO4_CO3",
            "Hydrotalcite_AlFe", "MSH",
        )
        @test n in names
    end
    @test all(length(end_members(p)) >= 2 for p in ss_all)
end
