# Chemical Kinetics

This tutorial covers the kinetics module of ChemistryLab.jl, which implements
mineral dissolution and precipitation kinetics coupled to fast aqueous speciation,
following the methodology of Leal et al. (2017).

## Background

The kinetics algorithm solves:

```math
\frac{d n_k}{dt} = -r_k(t), \quad k \in \text{kinetic minerals}
```

where ``r_k`` [mol/s] is the net dissolution rate of mineral ``k``.
At each evaluation of this ODE right-hand side, the aqueous speciation is
re-equilibrated, providing accurate activity
coefficients for the saturation ratio ``\Omega = \text{IAP}/K``.

The rate law follows transition-state theory (Palandri & Kharaka 2004):

```math
r_k = A_k \sum_m k_m(T) \cdot \prod_i a_i^{n_i} \cdot \text{sgn}(1-\Omega) \cdot |1-\Omega^p|^q
```

## Setting up rate constants

Rate constants are built as [`SymbolicFunc`](@ref) objects using the Arrhenius
factory, making them fully differentiable via ForwardDiff:

```julia
using ChemistryLab, ForwardDiff

# Arrhenius rate constant for calcite acid dissolution (Palandri & Kharaka 2004)
# Returns a NumericFunc — fully AD-compatible through k₀, Ea, and T
k_acid = arrhenius_rate_constant(5.012e-1, 14400.0)   # k₀ [mol/(m²s)], Ea [J/mol]

k_acid(; T = 298.15)    # → 0.5012  mol/(m²s)
k_acid(; T = 310.0)     # → higher value at elevated T

# AD through the rate constant
ForwardDiff.derivative(T -> k_acid(; T = T), 298.15)
```

For custom models, use [`add_kinetics_rate_model`](@ref) or build a
[`NumericFunc`](@ref) directly.

## Defining kinetic reactions

```julia
# Calcite acid + neutral mechanisms
k_acid    = arrhenius_rate_constant(5.012e-1, 14400.0)
k_neutral = arrhenius_rate_constant(1.549e-6, 23500.0)

model = TransitionStateRateModel([
    RateMechanism(k_acid,    1.0, 1.0, [RateModelCatalyst("H+", 1.0)]),
    RateMechanism(k_neutral, 1.0, 1.0),
])

# BET surface area: 0.09 m²/g × 1000 g/kg = 90 m²/kg
surface = BETSurfaceArea(90.0)

# idx_calcite: index of Calcite in the ChemicalSystem
# stoich_calcite: stoichiometric row vector
kr = KineticReaction(reaction_calcite, model, surface, idx_calcite, stoich_calcite)
```

## Running a simulation

```julia
using OrdinaryDiffEq   # activates KineticsOrdinaryDiffEqExt

# Build problem (0 to 2 hours)
kp = KineticsProblem(cs, [kr], state0, (0.0, 7200.0))

# Solve with default Rodas5P solver
sol = integrate(kp)

# Or with explicit solver and tolerances
ks  = KineticsSolver(; ode_solver=Rodas5P(), reltol=1e-8, abstol=1e-10)
sol = integrate(kp, ks)

# Access mineral moles at t = 1800 s
n_calcite_t = sol(1800.0)[1]
```

## Isothermal calorimetry

```julia
cal = IsothermalCalorimeter(298.15)   # T = 25 °C
sol = integrate(kp, ks; calorimeter = cal)

t, Q    = cumulative_heat(sol, cal)   # Q(t) [J]
t, qdot = heat_flow(sol, cal)         # q̇(t) [W]
```

## Semi-adiabatic calorimetry

The semi-adiabatic calorimeter solves

```math
\frac{dT}{dt} = \frac{\dot{q}(t) - \varphi(T(t) - T_{\rm env})}{C_p}
```

where `φ(ΔT)` is a user-supplied **heat-loss function**. The linear Newton's law
(`φ = L·ΔT`) is a convenient special case; the quadratic form used by
Lavergne et al. (2018) for cement hydration calorimetry is another common choice.

```julia
# Linear heat loss (Newton's law): backward-compatible L keyword
cal = SemiAdiabaticCalorimeter(;
    T0    = 293.15,   # initial temperature [K]
    T_env = 293.15,   # ambient temperature [K]
    L     = 0.5,      # heat-loss coefficient [W/K]
    Cp    = 4000.0,   # heat capacity of vessel + sample [J/K]
)

# Quadratic heat loss (Lavergne et al. 2018: φ(ΔT) = a·ΔT + b·ΔT²)
a, b = 0.48, 0.002   # [W/K] and [W/K²]
cal = SemiAdiabaticCalorimeter(;
    T0        = 293.15,
    T_env     = 293.15,
    Cp        = 4000.0,
    heat_loss = ΔT -> a * ΔT + b * ΔT^2,
)

sol = integrate(kp, ks; calorimeter = cal)

t, T_profile = temperature_profile(sol, cal)   # T(t) [K]
t, qdot      = heat_flow(sol, cal)             # q̇(t) [W]
t, Q         = cumulative_heat(sol, cal)       # Q(t) [J] (integrated from q̇)
```

## Parameter sensitivity and optimization

Because all rate constants are built as [`SymbolicFunc`](@ref)/[`NumericFunc`](@ref)
and the ODE function is ForwardDiff-compatible, you can differentiate the
solution trajectory with respect to kinetic parameters:

```julia
using ForwardDiff

# Sensitivity of final calcite amount w.r.t. activation energy
function n_calcite_final(Ea)
    k   = arrhenius_rate_constant(5.012e-1, Ea)
    rm  = TransitionStateRateModel([RateMechanism(k, 1.0, 1.0)])
    kr_ = KineticReaction(reaction_calcite, rm, surface, idx_calcite, stoich_calcite)
    kp_ = KineticsProblem(cs, [kr_], state0, (0.0, 3600.0))
    sol = integrate(kp_, ks)
    return sol.u[end][1]   # final moles of calcite
end

dndEa = ForwardDiff.derivative(n_calcite_final, 14400.0)
```

## Cement clinker hydration: Parrot–Killoh model

[`ParrotKillohRateModel`](@ref) implements the Parrot & Killoh (1984) kinetic model
for cement clinker hydration. Unlike surface-area-controlled models, it tracks the
degree of hydration α = 1 − n/n₀ and selects the controlling mechanism at each step
(nucleation-growth, interaction, or diffusion). The Arrhenius temperature correction
uses the activation energies of Schindler & Folliard (2005).

```julia
using ChemistryLab

# Predefined parameters for the four main clinker phases
rm_C3S  = PK_PARAMS_C3S    # alite       (K₁=1.5,  Ea=41 570 J/mol)
rm_C2S  = PK_PARAMS_C2S    # belite      (K₁=0.95, Ea=43 670 J/mol)
rm_C3A  = PK_PARAMS_C3A    # aluminate   (K₁=0.082,Ea=54 040 J/mol)
rm_C4AF = PK_PARAMS_C4AF   # ferrite     (K₁=0.165,Ea=34 420 J/mol)

# Evaluate at α = 0.3 (n_current = 0.7 mol, n_initial = 1.0 mol), T = 20 °C
r = PK_PARAMS_C3S(; T = 293.15, n_current = 0.7, n_initial = 1.0)  # mol/s

# Custom α_max (e.g. Powers 1948 limit for w/c = 0.40)
α_max = min(1.0, 0.40 / 0.42)
rm_custom = ParrotKillohRateModel(
    1.5, 3.3, 0.018, 2.5, 0.0024, 4.0, 0.5, 41_570.0;
    T_ref = 293.15, α_max = α_max,
)
```

The rate is AD-compatible:

```julia
using ForwardDiff

drdT = ForwardDiff.derivative(
    T -> PK_PARAMS_C3S(; T = T, n_current = 0.5, n_initial = 1.0),
    293.15,
)
```

## Full workflow: OPC clinker hydration with KineticsProblem

This example shows the complete chain — [`ChemicalSystem`](@ref) →
[`ChemicalState`](@ref) → [`KineticReaction`](@ref) → [`KineticsProblem`](@ref)
→ [`integrate`](@ref) — for a CEM I 52.5 R cement.

The convenience constructor `KineticReaction(cs, name, rate_model, surface_model)`
looks up the species by name and builds the stoichiometry automatically; no manual
index lookup is needed.

```julia
using ChemistryLab, OrdinaryDiffEq, DynamicQuantities, Printf

# ── 1. ChemicalSystem from the CEMDATA18 database ────────────────────────────
DATA_FILE = joinpath(pkgdir(ChemistryLab), "data", "cemdata18-thermofun.json")
substances = build_species(DATA_FILE)

input_species = split(
    "C3S C2S C3A C4AF " *
    "Portlandite Jennite ettringite monosulphate12 C3AH6 C3FH6 " *
    "H2O@",
)
species = speciation(substances, input_species; aggregate_state = [AS_AQUEOUS])
cs = ChemicalSystem(species, CEMDATA_PRIMARIES)

# ── 2. Initial state: 1 kg OPC (CEM I 52.5 R), w/c = 0.40 ───────────────────
const WC          = 0.40
const COMPOSITION = (C3S=0.619, C2S=0.165, C3A=0.080, C4AF=0.087)  # kg/kg cement

state0 = ChemicalState(cs)
for (name, frac) in pairs(COMPOSITION)
    set_quantity!(state0, string(name), frac * u"kg")
end
set_quantity!(state0, "H2O@", WC * u"kg")

# ── 3. Parrot–Killoh rate models with Powers α_max ───────────────────────────
α_max = min(1.0, WC / 0.42)

function pk_with_amax(template, α)
    ParrotKillohRateModel(
        template.K₁, template.N₁, template.K₂, template.N₂,
        template.K₃, template.N₃, template.B, template.Ea;
        T_ref = template.T_ref, α_max = α,
    )
end

pk_C3S  = pk_with_amax(PK_PARAMS_C3S,  α_max)
pk_C2S  = pk_with_amax(PK_PARAMS_C2S,  α_max)
pk_C3A  = pk_with_amax(PK_PARAMS_C3A,  α_max)
pk_C4AF = pk_with_amax(PK_PARAMS_C4AF, α_max)

# ── 4. Kinetic reactions: convenience constructor handles index + stoich ──────
# heat_per_mol [J/mol]: standard hydration enthalpies (Taylor 1997 / Lothenbach 2006)
# Required for SemiAdiabaticCalorimeter heat generation; ignored when calorimeter=nothing.
kr_C3S  = KineticReaction(cs, "C3S",  pk_C3S,  FixedSurfaceArea(1.0); heat_per_mol = 114_617.0)
kr_C2S  = KineticReaction(cs, "C2S",  pk_C2S,  FixedSurfaceArea(1.0); heat_per_mol =  44_782.0)
kr_C3A  = KineticReaction(cs, "C3A",  pk_C3A,  FixedSurfaceArea(1.0); heat_per_mol = 233_985.0)
kr_C4AF = KineticReaction(cs, "C4AF", pk_C4AF, FixedSurfaceArea(1.0); heat_per_mol = 203_617.0)

# ── 5. Problem + semi-adiabatic calorimeter ──────────────────────────────────
kp = KineticsProblem(
    cs, [kr_C3S, kr_C2S, kr_C3A, kr_C4AF], state0, (0.0, 7.0 * 86400.0);
    equilibrium_solver = nothing,
)

cal = SemiAdiabaticCalorimeter(;
    T0        = 293.15,
    T_env     = 293.15,
    Cp        = 1.0 * 800.0 + WC * 4186.0 + 1.0 * 900.0,   # ≈ 3449 J/K
    heat_loss = ΔT -> 0.30 * ΔT + 0.003 * ΔT^2,
)

# ── 6. Integrate and post-process ────────────────────────────────────────────
ks  = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1e-6, abstol = 1e-9)
sol = integrate(kp, ks; calorimeter = cal)

t_h = sol.t ./ 3600.0

_, T_vec = temperature_profile(sol, cal)
_, Q_vec = cumulative_heat(sol, cal)

n0 = sol.prob.p.n_initial
α_C3S  = 1.0 .- [u[1] for u in sol.u] ./ n0[1]
α_mean = (COMPOSITION.C3S .* α_C3S .+ ...) ./ sum(COMPOSITION)

@printf "ΔT_max = %.2f °C   Q_7d = %.1f kJ/kg\n" maximum(T_vec .- 273.15) - 20.0 Q_vec[end]/1000
```

!!! tip "Choosing `equilibrium_solver`"
    Setting `equilibrium_solver = nothing` skips the Gibbs minimization at each ODE
    step, which is appropriate for the Parrot–Killoh model (it does not use
    solution chemistry). For models that depend on the saturation ratio Ω — such as
    `TransitionStateRateModel` — pass an `EquilibriumSolver` to re-speciate the
    aqueous phase at every ODE evaluation.

!!! tip "Calorimetry with the convenience constructor"
    When `KineticReaction` is built via the convenience constructor (bare `Species`,
    no full reaction stoichiometry), the calorimeter cannot derive the heat of reaction
    from species formation enthalpies. Supply `heat_per_mol` [J/mol] explicitly — it
    takes priority over stoichiometric calculation and is the recommended approach for
    empirical models such as `ParrotKillohRateModel`. When `heat_per_mol = nothing`
    (default) and `kr.reaction isa AbstractReaction`, the enthalpy is computed
    from `Σ ν_j ΔₐH⁰_j` as usual.

## References

- Leal, A.M.M., Kulik, D.A., Smith, W.R., Saar, M.O. (2017).
  *An overview of computational methods for chemical equilibrium and kinetic
  calculations for geochemical and reactive transport modeling.*
  Pure and Applied Chemistry **89**, 597–643.
  <https://doi.org/10.1515/pac-2016-1107>

- Palandri, J.L. & Kharaka, Y.K. (2004).
  *A compilation of rate parameters of water-mineral interaction kinetics
  for application to geochemical modeling.*
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
