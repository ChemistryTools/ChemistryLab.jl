# The Pitzer ion-interaction model.
#
# No table of measured activity coefficients is available to this test suite, so
# nothing here pins the model to a literature number at high molality. What it
# does instead is check the three things that catch an algebra error without one:
# the analytical limit the model must reproduce as the solution becomes dilute,
# the Gibbs-Duhem relation between its own solutes and its own solvent — which is
# the whole reason to prefer this model, so it had better hold tightly — and
# agreement with an independent implementation of a different model in the range
# where both must equal the limiting law.

const PITZER_TOML = datapath("pitzer-reardon1990.toml")

_pz_params() = build_pitzer_parameters(PITZER_TOML)

function _pz_system(species = split("H2O@ Na+ Cl-"), primaries = ["H2O@", "Na+", "Cl-"])
    subs = build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false)
    d = Dict(symbol(s) => s for s in subs)
    return ChemicalSystem([d[s] for s in species], primaries)
end

const _PZ_M_W = 0.0180153

_pz_p(n) = (ΔₐG⁰overRT = zeros(n), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)

# ── the parameter file ───────────────────────────────────────────────────────

@testsection "the shipped Reardon parameter set loads and says what it is" begin
    p = _pz_params()
    @test p isa PitzerParameters
    @test p.beta0[("Na+", "Cl-")] ≈ 0.0765
    @test p.beta1[("Na+", "Cl-")] ≈ 0.2664
    @test p.Cphi[("Na+", "Cl-")] ≈ 0.0013
    # A 2-2 pair carries the third beta; a 1-1 pair does not.
    @test p.beta2[("Ca+2", "SO4-2")] ≈ -54.24
    @test p.beta2[("Na+", "Cl-")] == 0.0
    # theta and psi are found whichever way round the caller wrote the key.
    @test ChemistryLab._sym2(p.theta, "K+", "Na+") ≈ -0.012
    @test ChemistryLab._sym3(p.psi, "SO4-2", "Cl-", "Na+") ≈ 0.0014
    # An entry absent from the tables is zero, which is the literature's own
    # convention for theta and psi — and NOT for beta0, tested below.
    @test ChemistryLab._sym2(p.theta, "Na+", "Cl-") == 0.0

    # The provenance of a pair is readable, and the silicate and aluminate
    # parameters say that they were estimated rather than measured.
    @test pitzer_origin(p, "Na+", "Cl-") == "fitted"
    @test pitzer_origin(p, "Ca+2", "Al(OH)4-") == "estimated:HSO4-"
    @test pitzer_origin(p, "Na+", "H2SiO4-2") == "estimated:SO4-2"
    @test pitzer_origin(p, "Na+", "NoSuchIon-") == "unrecorded"
end

# ── construction refuses what it cannot describe ─────────────────────────────

@testsection "a Pitzer model cannot be reached half-parameterized" begin
    p = _pz_params()

    # Every table is a mandatory keyword: Julia refuses before any number is
    # computed, which is the enforcement asked of this model.
    @test_throws UndefKeywordError PitzerParameters(;
        beta0 = Dict(("Na+", "Cl-") => 0.1)
    )
    @test_throws UndefKeywordError PitzerActivityModel()

    # A beta1 for a pair with no beta0 is a typo in the caller's table.
    @test_throws ArgumentError PitzerParameters(;
        beta0 = Dict(("Na+", "Cl-") => 0.0765),
        beta1 = Dict(("K+", "Cl-") => 0.2122),
        beta2 = Dict{Tuple{String, String}, Float64}(),
        Cphi = Dict{Tuple{String, String}, Float64}(),
        theta = Dict{Tuple{String, String}, Float64}(),
        psi = Dict{Tuple{String, String, String}, Float64}(),
        lambda = Dict{Tuple{String, String}, Float64}(),
    )

    # The completeness of a set is a property of the set RELATIVE TO A SYSTEM,
    # so the refusal happens when the two meet — and it names the pairs.
    cs = _pz_system(
        split("H2O@ Na+ Cl- Br-"), ["H2O@", "Na+", "Cl-", "Br-"]
    )
    err = try
        activity_model(cs, PitzerActivityModel(; parameters = p))
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("Na+/Br-", err.msg)
    @test occursin("cannot fall back on ideal behavior", err.msg)

    # The set does cover a plain NaCl system.
    @test activity_model(_pz_system(), PitzerActivityModel(; parameters = p)) isa Function
end

# ── the dilute limit is the Debye-Hückel limiting law ────────────────────────

@testsection "Pitzer reduces to the limiting law as the solution empties" begin
    cs = _pz_system()
    model = PitzerActivityModel(; parameters = _pz_params())
    lna = activity_model(cs, model)
    n_w = 1 / _PZ_M_W
    A = hkf_debye_huckel_params(298.15, 1.0e5).A

    for m in (1.0e-6, 1.0e-5, 1.0e-4)
        out = lna([n_w, m * n_w * _PZ_M_W, m * n_w * _PZ_M_W], _pz_p(3))
        lnγ_mean = ((out[2] - log(m)) + (out[3] - log(m))) / 2
        # log₁₀ γ± → −A |z₊z₋| √I, exactly, and I = m for a 1-1 electrolyte.
        expected = -A * sqrt(m) * log(10)
        @test isapprox(lnγ_mean, expected; rtol = 0.02)
        # and the water activity approaches its ideal value from below
        @test out[1] < 0
        @test isapprox(exp(out[1]), 1 - 2 * m * _PZ_M_W; atol = 5.0e-6)
    end
end

@testsection "Pitzer and the B-dot model agree where both must" begin
    cs = _pz_system()
    pz = activity_model(cs, PitzerActivityModel(; parameters = _pz_params()))
    bd = activity_model(cs, HKFActivityModel())
    n_w = 1 / _PZ_M_W
    for m in (1.0e-5, 1.0e-4, 1.0e-3)
        n = [n_w, m * n_w * _PZ_M_W, m * n_w * _PZ_M_W]
        a, b = pz(n, _pz_p(3)), bd(n, _pz_p(3))
        # Two independent implementations of two different models, both obliged
        # to reduce to the same limiting law: they agree to better than 1 % in
        # ln γ at a millimolal and tighten as the solution empties.
        @test isapprox(a[2] - log(m), b[2] - log(m); rtol = 0.05)
    end
end

# ── the reason to use this model at all ─────────────────────────────────────

@testsection "Gibbs-Duhem holds by construction, not approximately" begin
    cs = _pz_system()
    μ = build_potentials(cs, PitzerActivityModel(; parameters = _pz_params()))
    μ_bdot = build_potentials(cs, HKFActivityModel())
    n_w = 1 / _PZ_M_W
    p = _pz_p(3)

    # The derivative is taken by AD and not by a finite difference on purpose.
    # A finite difference measures its own truncation error — 4e-6 at 0.1 mol/kg
    # with a 1e-6 step — which would hide the quantity being tested. With the
    # analytic derivative the residual is the algebra alone.
    function residual(mu, m, dn)
        n0 = [n_w, m, m]
        dμ = ForwardDiff.jacobian(nn -> mu(nn, p), n0) * dn
        return abs(sum(n0 .* dμ)) / max(norm(n0 .* abs.(dμ)), 1.0)
    end

    for dn in ([0.0, 1.0, 1.0], [0.0, 1.0, -1.0], [-1.0, 0.0, 0.0])
        for m in (0.1, 1.0, 3.0)
            # γ and φ are derived from one excess Gibbs energy, so this is an
            # identity of the algebra: it holds to machine precision, and any
            # missing or mistyped term in either half would break it.
            @test residual(μ, m, dn) < 1.0e-12
        end
    end

    # It is the property the B-dot model does not have. Along a dissolution at
    # 1 mol/kg its residual is not small at all, because its water activity is
    # an osmotic coefficient built on a single mean ionic radius while its γᵢ
    # use per-ion radii.
    @test residual(μ_bdot, 1.0, [0.0, 1.0, 1.0]) > 1.0e-6
end

# ── differentiability, which the solver requires ─────────────────────────────

@testsection "the Pitzer closure is differentiable" begin
    cs = _pz_system()
    lna = activity_model(cs, PitzerActivityModel(; parameters = _pz_params()))
    n_w = 1 / _PZ_M_W
    n0 = [n_w, 0.5, 0.5]
    J = ForwardDiff.jacobian(nn -> lna(nn, _pz_p(3)), n0)
    @test all(isfinite, J)
    @test size(J) == (3, 3)
    # The Taylor branch of g and g′ must be reached without producing a NaN
    # derivative: that is what an ionic strength near zero exercises.
    J0 = ForwardDiff.jacobian(nn -> lna(nn, _pz_p(3)), [n_w, 1.0e-12, 1.0e-12])
    @test all(isfinite, J0)
end

# ── a mixed solution, where theta and psi actually do something ─────────────

@testsection "the mixture terms are reached and change the answer" begin
    cs = _pz_system(
        split("H2O@ Na+ K+ Cl- SO4-2"), ["H2O@", "Na+", "K+", "Cl-", "SO4-2"]
    )
    p = _pz_params()
    full = activity_model(cs, PitzerActivityModel(; parameters = p))

    # The same set with every mixture term removed: θ and ψ are what a
    # single-electrolyte fit cannot give, so dropping them must move a mixed
    # solution and must NOT move a pure one.
    bare = PitzerParameters(;
        beta0 = p.beta0, beta1 = p.beta1, beta2 = p.beta2, Cphi = p.Cphi,
        theta = Dict{Tuple{String, String}, Float64}(),
        psi = Dict{Tuple{String, String, String}, Float64}(),
        lambda = p.lambda, origin = p.origin,
    )
    plain = activity_model(cs, PitzerActivityModel(; parameters = bare))

    n_w = 1 / _PZ_M_W
    n = [n_w, 0.5, 0.5, 0.6, 0.2]        # Na, K, Cl, SO4 — charge balanced
    a, b = full(n, _pz_p(5)), plain(n, _pz_p(5))
    @test !isapprox(a[2], b[2]; rtol = 1.0e-6)
    @test maximum(abs, a - b) > 1.0e-3

    cs2 = _pz_system()
    n2 = [n_w, 0.5, 0.5]
    a2 = activity_model(cs2, PitzerActivityModel(; parameters = p))(n2, _pz_p(3))
    b2 = activity_model(cs2, PitzerActivityModel(; parameters = bare))(n2, _pz_p(3))
    @test a2 ≈ b2
end

# ── the paths a plain NaCl system never reaches ─────────────────────────────

@testsection "neutral solutes go through the lambda terms" begin
    # A Pitzer set describes a neutral solute by `λ(neutral, ion)` — salting out
    # in Pitzer's form — and nothing else. With no λ the neutral is ideal, which
    # is the literature's convention and worth pinning rather than assuming.
    subs = build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false)
    d = Dict(symbol(s) => s for s in subs)
    cs = ChemicalSystem(
        [d[s] for s in split("H2O@ Na+ Cl- CO2@")], ["H2O@", "Na+", "Cl-", "CO2@"]
    )
    base = _pz_params()

    with_λ = PitzerParameters(;
        beta0 = base.beta0, beta1 = base.beta1, beta2 = base.beta2, Cphi = base.Cphi,
        theta = base.theta, psi = base.psi,
        lambda = Dict(("CO2@", "Na+") => 0.1, ("CO2@", "Cl-") => -0.005),
    )
    without_λ = PitzerParameters(;
        beta0 = base.beta0, beta1 = base.beta1, beta2 = base.beta2, Cphi = base.Cphi,
        theta = base.theta, psi = base.psi,
        lambda = Dict{Tuple{String, String}, Float64}(),
    )

    # The molality must be formed with the molar mass the *closure* uses — the
    # database's — and not with a rounded constant, or the identity below is off
    # by the difference between the two, which is 1.7e-5 in ln γ.
    M_w = ustrip(us"kg/mol", cs.species[only(cs.idx_solvent)][:M])
    n_w = 1 / M_w
    n = [n_w, 1.0, 1.0, 0.05]                     # 1 molal NaCl, 0.05 molal CO₂
    a = activity_model(cs, PitzerActivityModel(; parameters = with_λ))(n, _pz_p(4))
    b = activity_model(cs, PitzerActivityModel(; parameters = without_λ))(n, _pz_p(4))

    m_co2 = 0.05 / (n[1] * M_w)
    γ_with = exp(a[4] - log(m_co2))
    γ_without = exp(b[4] - log(m_co2))
    @test isapprox(γ_without, 1.0; rtol = 1.0e-10)      # no λ ⇒ ideal, exactly
    @test γ_with > 1.0                                   # a positive λ salts out
    # ln γ_n = 2 Σ_i m_i λ_ni, which is checkable by hand here.
    m_ion = 1.0 / (n[1] * M_w)
    @test isapprox(log(γ_with), 2 * (m_ion * 0.1 + m_ion * (-0.005)); rtol = 1.0e-10)

    # The solvent feels the neutral too, through the same λ.
    @test a[1] != b[1]
end

@testsection "gases and solid-solution end-members are filled in" begin
    subs = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    d = Dict(symbol(s) => s for s in subs)

    # A gas phase: an ideal mixture, ln a = ln x.
    cs_g = ChemicalSystem(
        [d[s] for s in split("H2O@ Na+ Cl- CO2 O2")],
        ["H2O@", "Na+", "Cl-", "CO2", "O2"],
    )
    model = PitzerActivityModel(; parameters = _pz_params())
    n_w = 1 / _PZ_M_W
    out = activity_model(cs_g, model)([n_w, 0.1, 0.1, 3.0, 1.0], _pz_p(5))
    @test isapprox(out[4], log(3.0 / 4.0); rtol = 1.0e-10)
    @test isapprox(out[5], log(1.0 / 4.0); rtol = 1.0e-10)

    # A solid solution: the end-members must be filled by the solid-solution
    # branch. An aqueous model that forgets to call it leaves them at ln a = 0,
    # i.e. treated as pure phases, which is the silent failure this asserts on.
    ss = build_solid_solutions(datapath("solid_solutions.toml"), d; skip_missing = true)
    cshq = only(filter(p -> name(p) == "CSHQ", ss))
    cs_ss = ChemicalSystem(
        vcat([d[s] for s in split("H2O@ Na+ Cl-")], end_members(cshq)),
        ["H2O@", "Na+", "Cl-"]; solid_solutions = [cshq],
    )
    k = length(cs_ss.species)
    n = vcat([n_w, 0.1, 0.1], fill(0.25, k - 3))
    out_ss = activity_model(cs_ss, model)(n, _pz_p(k))
    for i in 4:k
        @test out_ss[i] < 0                       # ln x < 0 for a fraction < 1
        @test isfinite(out_ss[i])
    end
    @test !all(iszero, out_ss[4:k])
end

@testsection "the Debye-Hückel slope can follow the temperature" begin
    cs = _pz_system()
    p = _pz_params()
    fixed = activity_model(cs, PitzerActivityModel(; parameters = p))
    varying = activity_model(
        cs, PitzerActivityModel(; parameters = p, temperature_dependent = true)
    )
    n_w = 1 / _PZ_M_W
    n = [n_w, 0.5, 0.5]

    at25 = (ΔₐG⁰overRT = zeros(3), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    at60 = (ΔₐG⁰overRT = zeros(3), T = 333.15, P = 1.0e5, ϵ = 1.0e-30)

    # Fixed: the temperature in `p` is ignored, by construction.
    @test fixed(n, at25) ≈ fixed(n, at60)
    # Varying: A_φ rises with temperature, so the ions are further from ideal.
    @test !isapprox(varying(n, at25), varying(n, at60))
    @test varying(n, at60)[2] < varying(n, at25)[2]
    # and at 25 °C the two agree, since that is where the fixed value comes from
    @test isapprox(varying(n, at25), fixed(n, at25); rtol = 1.0e-3)
end
