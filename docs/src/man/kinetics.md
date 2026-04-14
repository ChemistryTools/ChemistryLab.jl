# Chemical Kinetics

This tutorial covers the kinetics module of ChemistryLab.jl, which implements
mineral dissolution and precipitation kinetics coupled to fast aqueous speciation,
following the methodology of Leal et al. (2017).

## Background

The kinetics algorithm solves:

```math
\frac{d n_k}{dt} = r_k(t), \quad k \in \text{kinetic minerals}
```

where ``r_k`` [mol/s] is the net rate of reaction ``k`` (positive = dissolution,
negative stoichiometric coefficient for the mineral so `dn/dt < 0`).
At each evaluation of this ODE right-hand side, the aqueous speciation is optionally
re-equilibrated, providing accurate activity coefficients for the saturation ratio
``\Omega = \text{IAP}/K``.

## Rate functions: KineticFunc and StateView

All rate functions are encapsulated as [`KineticFunc`](@ref) objects.  They share the
**positional six-argument calling convention**:

```julia
kf(T, P, t, n::StateView, lna::StateView, n_initial::StateView) -> Real [mol/s]
```

where

| Argument | Unit | Description |
| --- | --- | --- |
| `T` | K | Temperature (plain `Real` or `ForwardDiff.Dual`) |
| `P` | Pa | Pressure |
| `t` | s | Current integration time |
| `n` | mol | Moles of all species — named access via [`StateView`](@ref) |
| `lna` | — | Log-activities — named access via `StateView` |
| `n_initial` | mol | Initial moles — named access via `StateView` (always `Float64`) |

[`StateView`](@ref) provides O(1) named access to a species vector via a pre-built
dictionary (`sv["C3S"]` — no per-step dict allocation):

```julia
index = Dict("C3S" => 1, "C2S" => 2)
data  = [2.5, 1.8]          # moles
sv    = StateView(data, index)

sv["C3S"]          # → 2.5
haskey(sv, "C3S")  # → true
```

The `index` dict is built **once** at `KineticsProblem` construction; the same dict is
shared by all `StateView` wrappers inside the ODE hot path.

## Setting up rate constants

Rate constants are built as [`NumericFunc`](@ref) objects using the Arrhenius factory:

```julia
using ChemistryLab, ForwardDiff

# Arrhenius rate constant for calcite acid dissolution (Palandri & Kharaka 2004)
k_acid = arrhenius_rate_constant(5.012e-1, 14400.0)   # k₀ [mol/(m²s)], Ea [J/mol]

k_acid(; T = 298.15)    # → 0.5012 mol/(m²s)
k_acid(; T = 310.0)     # → higher value at elevated T

ForwardDiff.derivative(T -> k_acid(; T = T), 298.15)  # AD-compatible
```

## Cement clinker hydration: Parrot–Killoh model

[`parrot_killoh`](@ref) is the factory for the Parrot & Killoh (1984) kinetic model.
It returns a [`KineticFunc`](@ref) that uses the `StateView`-based calling convention
and is fully AD-compatible.

```julia
using ChemistryLab, DynamicQuantities

# Predefined NamedTuple parameters for the four main clinker phases
#   PK_PARAMS_C3S, PK_PARAMS_C2S, PK_PARAMS_C3A, PK_PARAMS_C4AF
# Each has keys: K₁, N₁, K₂, N₂, K₃, N₃, B, Ea, T_ref (with units)

pk_C3S = parrot_killoh(PK_PARAMS_C3S, "C3S")    # → KineticFunc

# Custom α_max (Powers 1948 limit for w/c = 0.40)
WC    = 0.40
α_max = min(1.0, WC / 0.42)
pk_C3S_wc = parrot_killoh(PK_PARAMS_C3S, "C3S"; α_max = α_max)

# Custom parameters — accepts Quantity or plain SI Real
pk_ggbs = parrot_killoh(
    (K₁=0.15u"1/d", N₁=2.0, K₂=0.003u"1/d", N₂=2.0,
     K₃=0.0015u"1/d", N₃=3.5, B=0.2, Ea=46_000.0u"J/mol", T_ref=293.15u"K"),
    "GGBS";
    α_max = 0.9,
)
```

Rate evaluation:

```julia
index = Dict("C3S" => 1)
n_sv  = StateView([0.5], index)    # 0.5 mol C3S remaining
n0_sv = StateView([1.0], index)    # 1.0 mol initial
lna   = StateView([0.0], index)

r = pk_C3S(293.15, 1.0e5, 0.0, n_sv, lna, n0_sv)   # mol/s
```

AD smoke-test:

```julia
using ForwardDiff

drdT = ForwardDiff.derivative(T -> pk_C3S(T, 1.0e5, 0.0, n_sv, lna, n0_sv), 293.15)
@assert isfinite(drdT) && drdT > 0
```

## Transition-state theory rate model

For models based on solution chemistry (calcite, quartz, …),
[`transition_state`](@ref) builds a multi-mechanism TST [`KineticFunc`](@ref):

```julia
k_neutral = arrhenius_rate_constant(1.549e-6, 23500.0)
k_acid    = arrhenius_rate_constant(5.012e-1, 14400.0)

surface = BETSurfaceArea(90.0)   # 90 m²/kg

tst = transition_state(
    [
        RateMechanism(k_neutral, 1.0, 1.0),
        RateMechanism(k_acid,    1.0, 1.0, [RateModelCatalyst("H+", 1.0)]),
    ],
    cs,
    cs.dict_reactions["calcite dissolution"],
    surface,
)

# tst is a KineticFunc — same calling convention as parrot_killoh
kr = KineticReaction(cs, cs.dict_reactions["calcite dissolution"], tst)
```

The `transition_state` factory captures the ΔₐG⁰ callables from the reaction and
recomputes `ln K(T)` at every ODE step, so the saturation ratio Ω is correct even
when the temperature changes (semi-adiabatic calorimetry).

For a single mechanism without catalysts, use [`first_order_rate`](@ref) as a
convenience wrapper.

## Defining kinetic reactions

[`KineticReaction`](@ref) associates a [`Reaction`](@ref) with a `KineticFunc`:

```julia
pk = parrot_killoh(PK_PARAMS_C3S, "C3S")

# Convenience: look up "C3S" in cs by symbol/formula, build minimal dissolution Reaction
kr = KineticReaction(cs, "C3S", pk)

# Explicit Reaction object (e.g. from cs.dict_reactions)
rxn = cs.dict_reactions["calcite dissolution"]
kr  = KineticReaction(cs, rxn, tst_calcite)

# Reaction-centric (kinetics stored in rxn.properties)
rxn[:rate] = parrot_killoh(PK_PARAMS_C3S, "C3S")
kr = KineticReaction(cs, rxn)
```

## [Reaction-centric workflow](@id reaction-centric-workflow)

Kinetics can be attached directly to `Reaction` objects via their `properties` dict:

### Properties table

| Key | Type | Required | Description |
| --- | --- | :---: | --- |
| `:rate` | `KineticFunc` or `(T,P,t,n,lna,n₀)->Real` | ✓ | Net dissolution rate [mol/s] |
| `:heat_per_mol` | `Number` | | Override molar heat [J/mol], positive = exothermic; needed only for custom species without `ΔₐH⁰` |

When `:rate` is a plain callable (not a `KineticFunc`), it is automatically wrapped
in a `KineticFunc` with empty `refs`.

### Example: C₃S hydration (reaction-centric)

```julia
using ChemistryLab, OrdinaryDiffEq, DynamicQuantities

sp(name) = cs[name]

pk_C3S = parrot_killoh(PK_PARAMS_C3S, "C3S"; α_max = min(1.0, WC / 0.42))

rxn_C3S = Reaction(
    OrderedDict(sp("C3S") => 1.0, sp("H2O@") => 3.33),
    OrderedDict(sp("Jennite") => 0.167, sp("Portlandite") => 1.5);
    symbol = "C₃S hydration",
)
rxn_C3S[:rate] = pk_C3S   # heat computed from species ΔₐH⁰ (CEMDATA18)

kp = KineticsProblem(cs, [rxn_C3S, rxn_C2S, rxn_C3A, rxn_C4AF], state0, tspan)
```

## Unit-aware API

All constructors support DynamicQuantities `Quantity` values.  Plain `Real` → SI assumed.

| Constructor | Parameter | Default SI unit | Example |
| --- | --- | --- | --- |
| `arrhenius_rate_constant` | `k₀` | mol/(m²·s) | `0.5u"mmol/(m^2*s)"` |
| `arrhenius_rate_constant` | `Ea` | J/mol | `62.0u"kJ/mol"` |
| `parrot_killoh` params | `K₁, K₂, K₃` | s⁻¹ | `1.5u"1/d"` |
| `parrot_killoh` params | `Ea` | J/mol | `41.57u"kJ/mol"` |
| `parrot_killoh` params | `T_ref` | K | `293.15u"K"` |
| `FixedSurfaceArea` | `A` | m² | `500.0u"cm^2"` |
| `BETSurfaceArea` | `A_specific` | m²/kg | `0.09u"m^2/g"` |
| `KineticReaction` | `heat_per_mol` (custom species) | J/mol | `36.1u"kJ/mol"` |
| `KineticsProblem` | `tspan` | s | `(0.0u"s", 7.0u"d")` |
| `IsothermalCalorimeter` | `T` | K | `298.15u"K"` |
| `SemiAdiabaticCalorimeter` | `Cp` | J/K | `4.0u"kJ/K"` |
| `SemiAdiabaticCalorimeter` | `T_env` | K | `293.15u"K"` |
| `SemiAdiabaticCalorimeter` | `L` | W/K | `500.0u"mW/K"` |
| `SemiAdiabaticCalorimeter` | `T0` | K | `293.15u"K"` |

## Running a simulation

```julia
using OrdinaryDiffEq   # activates KineticsOrdinaryDiffEqExt

kp  = KineticsProblem(cs, [kr], state0, (0.0, 7200.0))
sol = integrate(kp)   # default Rodas5P

ks  = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1e-8, abstol = 1e-10)
sol = integrate(kp, ks)
```

## Isothermal calorimetry

```julia
cal = IsothermalCalorimeter(298.15u"K")   # T = 25 °C
kp = KineticsProblem(cs, reactions, state0, tspan; calorimeter = cal)
sol = integrate(kp, ks)

t, Q    = cumulative_heat(sol, cal)   # Q(t) [J]
t, qdot = heat_flow(sol, cal)         # q̇(t) [W]
```

## Semi-adiabatic calorimetry (Lavergne et al. 2018)

The semi-adiabatic calorimeter solves:

```math
\frac{dT}{dt} = \frac{\dot{q}(t) - \varphi(T(t) - T_{\rm env})}{C_p + \sum_i n_i C^\circ_{p,i}(T)}
```

The denominator uses the **variable total heat capacity** `Cp_total = Cp + Σᵢ nᵢ Cp°ᵢ(T)`,
where `Cp°ᵢ(T)` are the molar heat capacities from the thermodynamic database
(Lavergne et al. 2018).

[`SemiAdiabaticCalorimeter`](@ref) bundles hardware parameters and initial temperature:

```julia
using ChemistryLab, DynamicQuantities

# Quadratic heat loss (Lavergne et al. 2018: φ = a·ΔT + b·ΔT²)
cal = SemiAdiabaticCalorimeter(;
    Cp        = (1.0 * 800.0 + WC * 4186.0 + 1.0 * 900.0) * u"J/K",   # ≈ 3449 J/K
    T_env     = 293.15u"K",
    heat_loss = ΔT -> 0.3 * ΔT + 0.003 * ΔT^2,
    T0        = 293.15u"K",
)

# Shorthand for linear Newton cooling (L keyword)
cal_lin = SemiAdiabaticCalorimeter(;
    Cp = 4000.0u"J/K", T_env = 293.15u"K", L = 0.5u"W/K", T0 = 293.15u"K",
)

kp = KineticsProblem(cs, reactions, state0, tspan; calorimeter = cal)
sol = integrate(kp, ks)

t, T_profile = temperature_profile(sol, cal)   # T(t) [K]
t, qdot      = heat_flow(sol, cal)             # q̇(t) [W]
t, Q         = cumulative_heat(sol, cal)       # Q(t) [J]
```

## Full workflow: OPC clinker hydration with KineticsProblem

Complete chain — [`ChemicalSystem`](@ref) → [`ChemicalState`](@ref) →
[`KineticReaction`](@ref) → [`KineticsProblem`](@ref) → [`integrate`](@ref):

```julia
using ChemistryLab, OrdinaryDiffEq, DynamicQuantities, Printf

# ── 1. ChemicalSystem from CEMDATA18 ────────────────────────────────────────
DATA_FILE = joinpath(pkgdir(ChemistryLab), "data", "cemdata18-thermofun.json")
substances = build_species(DATA_FILE)

input_species = split(
    "C3S C2S C3A C4AF " *
    "Portlandite Jennite ettringite monosulphate12 C3AH6 C3FH6 H2O@",
)
species = speciation(substances, input_species; aggregate_state = [AS_AQUEOUS])
cs = ChemicalSystem(species, CEMDATA_PRIMARIES)

# ── 2. Initial state: 1 kg OPC (CEM I 52.5 R), w/c = 0.40 ─────────────────
WC          = 0.40
COMPOSITION = (C3S=0.619, C2S=0.165, C3A=0.080, C4AF=0.087)

state0 = ChemicalState(cs)
for (name, frac) in pairs(COMPOSITION)
    set_quantity!(state0, string(name), frac * u"kg")
end
set_quantity!(state0, "H2O@", WC * u"kg")

# ── 3. Parrot–Killoh rate functions with Powers α_max ───────────────────────
α_max   = min(1.0, WC / 0.42)
pk_C3S  = parrot_killoh(PK_PARAMS_C3S,  "C3S";  α_max)
pk_C2S  = parrot_killoh(PK_PARAMS_C2S,  "C2S";  α_max)
pk_C3A  = parrot_killoh(PK_PARAMS_C3A,  "C3A";  α_max)
pk_C4AF = parrot_killoh(PK_PARAMS_C4AF, "C4AF"; α_max)

# ── 4. Kinetic reactions (reaction-centric) ─────────────────────────────────
# Reactions follow Lothenbach & Winnefeld (2006) — Jennite = Ca₉Si₆O₁₈(OH)₆·8H₂O
# Calorimetric heat is computed automatically from CEMDATA18 species ΔₐH⁰.
sp(name) = cs[name]

rxn_C3S = Reaction(
    OrderedDict(sp("C3S") => 1.0, sp("H2O@") => 3.33),
    OrderedDict(sp("Jennite") => 0.167, sp("Portlandite") => 1.5);
    symbol = "C₃S hydration",
)
rxn_C3S[:rate] = pk_C3S

rxn_C2S = Reaction(
    OrderedDict(sp("C2S") => 1.0, sp("H2O@") => 2.33),
    OrderedDict(sp("Jennite") => 0.167, sp("Portlandite") => 0.5);
    symbol = "C₂S hydration",
)
rxn_C2S[:rate] = pk_C2S

rxn_C3A = Reaction(
    OrderedDict(sp("C3A") => 1.0, sp("H2O@") => 6.0),
    OrderedDict(sp("C3AH6") => 1.0);
    symbol = "C₃A hydration",
)
rxn_C3A[:rate] = pk_C3A

rxn_C4AF = Reaction(
    OrderedDict(sp("C4AF") => 1.0, sp("H2O@") => 6.0),
    OrderedDict(sp("C3AH6") => 0.5, sp("C3FH6") => 0.5, sp("Portlandite") => 1.0);
    symbol = "C₄AF hydration",
)
rxn_C4AF[:rate] = pk_C4AF

# ── 5. Problem + semi-adiabatic calorimeter ─────────────────────────────────
cal = SemiAdiabaticCalorimeter(;
    Cp        = (1.0 * 800.0 + WC * 4186.0 + 1.0 * 900.0) * u"J/K",
    T_env     = 293.15u"K",
    heat_loss = ΔT -> 0.30 * ΔT + 0.003 * ΔT^2,
    T0        = 293.15u"K",
)

kp = KineticsProblem(
    cs, [rxn_C3S, rxn_C2S, rxn_C3A, rxn_C4AF], state0, (0.0, 7.0 * 86400.0);
    calorimeter = cal,
    equilibrium_solver = nothing,
)

# ── 6. Integrate and post-process ───────────────────────────────────────────
ks  = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1e-6, abstol = 1e-9)
sol = integrate(kp, ks)

_, T_vec = temperature_profile(sol, cal)
_, Q_vec = cumulative_heat(sol, cal)

n0_kin = [sol.prob.p.n_initial_full[i] for i in kp.idx_kinetic]
n_kin  = [[u[i] for u in sol.u] for i in eachindex(n0_kin)]

function phase_alpha(cs, kp, n0_kin, n_kin, name)
    sp_idx = findfirst(s -> ChemistryLab.phreeqc(ChemistryLab.formula(s)) == name, cs.species)
    pos    = findfirst(==(sp_idx), kp.idx_kinetic)
    isnothing(pos) && return fill(NaN, length(sol.t))
    return 1.0 .- n_kin[pos] ./ n0_kin[pos]
end

α_C3S = phase_alpha(cs, kp, n0_kin, n_kin, "C3S")

@printf "ΔT_max = %.2f °C   Q_7d = %.1f kJ/kg   α(C3S) = %.4f\n" \
    maximum(T_vec .- 273.15) - 20.0  Q_vec[end]/1000  α_C3S[end]
```

!!! tip "Choosing `equilibrium_solver`"
    Setting `equilibrium_solver = nothing` skips the Gibbs minimization at each ODE
    step, which is appropriate for the Parrot–Killoh model (it does not use
    solution chemistry).  For `transition_state` models, pass an `EquilibriumSolver`
    to re-speciate the aqueous phase at every ODE evaluation.

!!! tip "Calorimetry and ΔᵣH⁰"
    The calorimeter computes the heat generation rate as `q̇ = Σ rᵢ × (−ΔᵣH⁰ᵢ)`,
    where `ΔᵣH⁰ᵢ(T)` is the reaction enthalpy built automatically from species
    `ΔₐH⁰` properties via `complete_thermo_functions!`.  For database species
    (e.g. CEMDATA18) no extra input is needed.  For **custom species** that lack
    a `ΔₐH⁰` entry (GGBS, MK, …), set `:heat_per_mol` [J/mol] on the reaction
    to override (positive = exothermic).

## Multiple kinetic reactions per mineral

Following Leal et al. (2017), **reactions** — not species — carry kinetics.
A single mineral can appear in multiple `KineticReaction` objects
(e.g. C₃A via an early ettringite pathway and a late monosulphate pathway).
The ODE state contains **one entry per unique mineral**; contributions accumulate
in `du[j]`.

```julia
pk_c3a = parrot_killoh(PK_PARAMS_C3A, "C3A"; α_max)

kr_C3A_ett  = KineticReaction(cs, rxn_C3A_ettringite,   pk_c3a)
kr_C3A_mono = KineticReaction(cs, rxn_C3A_monosulphate, pk_c3a)

kp = KineticsProblem(
    cs, [kr_C3S, kr_C2S, kr_C3A_ett, kr_C3A_mono, kr_C4AF], state0, tspan,
)
# length(build_u0(kp)) == 4  (one entry per unique mineral: C3S, C2S, C3A, C4AF)
```

!!! note "Species classification (Leal et al. 2015)"
    `KineticsProblem` automatically applies the Leal et al. 2015 classification:
    - **Kinetic species** (tracked in ODE state `u`): `AS_CRYSTAL` species with non-zero
      stoichiometry in any kinetic reaction.
    - **Equilibrium species**: aqueous species re-equilibrated by `equilibrium_solver`.
    - **Inert species**: all others.

## Parameter sensitivity and optimization

The entire chain is ForwardDiff-compatible:

```julia
using ForwardDiff

function n_C3S_final(K₁)
    params = merge(PK_PARAMS_C3S, (K₁ = K₁,))   # override K₁ (plain Real, SI = 1/s)
    pk  = parrot_killoh(params, "C3S")
    kr_ = KineticReaction(cs, "C3S", pk)
    kp_ = KineticsProblem(cs, [kr_], state0, (0.0, 86400.0))
    sol = integrate(kp_, ks)
    return sol.u[end][1]
end

dn_dK₁ = ForwardDiff.derivative(n_C3S_final, safe_ustrip(us"1/s", PK_PARAMS_C3S.K₁))
```

## References

- Leal, A.M.M., Kulik, D.A., Smith, W.R., Saar, M.O. (2017).
  *An overview of computational methods for chemical equilibrium and kinetic
  calculations for geochemical and reactive transport modeling.*
  Pure and Applied Chemistry **89**, 597–643.
  <https://doi.org/10.1515/pac-2016-1107>

- Palandri, J.L. & Kharaka, Y.K. (2004).
  *A compilation of rate parameters of water-mineral interaction kinetics.*
  USGS Open-File Report 2004-1068.

- Parrot, L.J. & Killoh, D.C. (1984).
  *Prediction of cement hydration.*
  British Ceramic Proceedings **35**, 41–53.

- Lavergne, F., Ben Fraj, A., Bayane, I., Barthélémy, J.-F. (2018).
  *Estimating the mechanical properties of hydrating blended cementitious materials.*
  Cement and Concrete Research **104**, 37–60.
  <https://doi.org/10.1016/j.cemconres.2017.11.007>

- Schindler, A.K. & Folliard, K.J. (2005).
  *Heat of hydration models for cementitious materials.*
  ACI Materials Journal **102**, 24–33.
