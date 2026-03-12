# [Solid Solutions](@id sec-solid-solutions)

A **solid solution** is a crystalline phase in which several **end-members** mix
at the atomic scale.  Classic examples in geochemistry and cement chemistry include:

| Solid solution | End-members |
|:---|:---|
| Calcite–Rhodochrosite | CaCO₃ + MnCO₃ |
| Calcite–Magnesite | CaCO₃ + MgCO₃ |
| Hydrocalumite (AFm-SO₄/CO₃) | Ca₄Al₂(OH)₁₂·SO₄ + Ca₄Al₂(OH)₁₂·CO₃ |
| C-S-H | variable Ca/Si ratio |

Unlike a **pure mineral** (activity = 1), the activity of end-member $i$ in a solid
solution depends on its mole fraction $x_i = n_i / \sum_j n_j$:

$$\ln a_i = \ln x_i + \ln \gamma_i$$

where $\ln \gamma_i = 0$ for ideal mixing and is given by a Margules expression for
non-ideal (regular) mixing.

---

## Mixing models

Two models are available, both implementing the
[`AbstractSolidSolutionMixingModel`](@ref) interface.

### Ideal mixing — `IdealSolidSolutionMixingModel`

$$\ln a_i = \ln x_i, \qquad \ln \gamma_i = 0$$

### Symmetric Margules — `RegularSolidSolutionMixingModel`

The excess Gibbs energy per mole of phase is:

$$G^{xs} = \sum_{i < j} W_{ij}\, x_i\, x_j$$

where $W_{ij}$ [J/mol] is the symmetric Margules interaction parameter between
end-members $i$ and $j$.  The activity coefficients are the exact derivatives of
$n\,G^{xs}$ with respect to $n_i$:

$$RT \ln \gamma_i =
    \sum_{j \neq i} W_{ij}\, x_j (1 - x_i)
  - \sum_{\substack{j < k \\ j \neq i,\, k \neq i}} W_{jk}\, x_j\, x_k$$

For a **binary** system this simplifies to
$RT \ln \gamma_1 = W_{12}\, x_2^2$ and $RT \ln \gamma_2 = W_{12}\, x_1^2$.

---

## Example: Ca–Mg carbonate solid solution

We study the dissolution of a calcite–magnesite solid solution in water:

$$(\text{Ca}_{1-x}\text{Mg}_x)\text{CO}_3 \;\rightleftharpoons\;
  (1-x)\,\text{Ca}^{2+} + x\,\text{Mg}^{2+} + \text{CO}_3^{2-}$$

The two end-members are **calcite** (`Cal`, CaCO₃) and **magnesite** (`Mgs`, MgCO₃).

### Loading species from the SLOP98 database

All species are loaded from the **SLOP98 inorganic** database.

| Symbol   | Species                | Phase          |
|:---------|:-----------------------|:---------------|
| `H2O@`   | H₂O — water            | aqueous solvent |
| `H+`     | H⁺ — proton            | aqueous solute  |
| `OH-`    | OH⁻ — hydroxide        | aqueous solute  |
| `Ca+2`   | Ca²⁺                   | aqueous solute  |
| `Mg+2`   | Mg²⁺                   | aqueous solute  |
| `CO3-2`  | CO₃²⁻ — carbonate      | aqueous solute  |
| `HCO3-`  | HCO₃⁻ — bicarbonate    | aqueous solute  |
| `Cal`    | CaCO₃ — calcite        | crystal         |
| `Mgs`    | MgCO₃ — magnesite      | crystal         |

```@setup ss_example
using ChemistryLab
using DynamicQuantities

substances = build_species("../../../data/slop98-inorganic-thermofun.json")
dict = Dict(symbol(s) => s for s in substances)

species = [dict[sym] for sym in split("H2O@ H+ OH- Ca+2 Mg+2 CO3-2 HCO3- Cal Mgs")]
cs = ChemicalSystem(species, ["H2O@", "H+", "CO3-2", "Ca+2", "Mg+2"])
```

```julia
using ChemistryLab
using DynamicQuantities

substances = build_species("../../../data/slop98-inorganic-thermofun.json")
dict = Dict(symbol(s) => s for s in substances)

species = [dict[sym] for sym in split("H2O@ H+ OH- Ca+2 Mg+2 CO3-2 HCO3- Cal Mgs")]
cs = ChemicalSystem(species, ["H2O@", "H+", "CO3-2", "Ca+2", "Mg+2"])
```

### Reference solubility products

The dissolution reactions for each pure end-member are:

$$\text{Cal (CaCO}_3\text{)} \rightleftharpoons \text{Ca}^{2+} + \text{CO}_3^{2-}
  \qquad \log K_{\text{sp,Cal}} \approx -8.48$$

$$\text{Mgs (MgCO}_3\text{)} \rightleftharpoons \text{Mg}^{2+} + \text{CO}_3^{2-}
  \qquad \log K_{\text{sp,Mgs}} \approx -8.03$$

```@example ss_example
T = 298.15u"K"
R = ustrip(Constants.R)  # J/(mol·K)
RT = R * ustrip(us"K", T)

sp = Dict(symbol(s) => s for s in cs.species)

# Calcite dissolution: Cal → Ca²⁺ + CO₃²⁻
ΔrG_cal = sp["Ca+2"].ΔₐG⁰(T=T) + sp["CO3-2"].ΔₐG⁰(T=T) - sp["Cal"].ΔₐG⁰(T=T)
logKsp_cal = -ΔrG_cal / (RT * log(10))

# Magnesite dissolution: Mgs → Mg²⁺ + CO₃²⁻
ΔrG_mag = sp["Mg+2"].ΔₐG⁰(T=T) + sp["CO3-2"].ΔₐG⁰(T=T) - sp["Mgs"].ΔₐG⁰(T=T)
logKsp_mag = -ΔrG_mag / (RT * log(10))

@show logKsp_cal   # ≈ −8.48  (SLOP98)
@show logKsp_mag   # ≈ −8.03  (SLOP98)
nothing
```

### Effect of ideal mixing on the effective solubility product

When both minerals form an ideal solid solution, the effective dissolution equilibrium for
end-member $i$ becomes:

$$\log K_{\text{eff},i} = \log K_{\text{sp},i} - \log_{10}(x_i)$$

Because $x_i \leq 1$, the solid solution is **less soluble** than the pure end-member.

```@example ss_example
using Printf

println("x_Cal  |  log Keff (Cal)  |  x_Mgs   |  log Keff (Mgs)")
println("-------|------------------|----------------")
for x_cal in [0.9, 0.7, 0.5, 0.3, 0.1]
    x_mag = 1.0 - x_cal
    @printf("  %.1f  |     %+.2f        |     %.1f  |     %+.2f\n",
            x_cal,
            logKsp_cal - log10(max(x_cal, 1e-10)),
            x_mag,
            logKsp_mag - log10(max(x_mag, 1e-10)))
end
```

---

### Using `SolidSolutionActivityModel` in an equilibrium calculation

A solid solution is introduced by wrapping any base activity model with
[`SolidSolutionActivityModel`](@ref).
The end-member symbols passed to [`SolidSolution`](@ref) must match the
`symbol` of the species in the `ChemicalSystem` (here `"Cal"` and `"Mgs"`).

```@example ss_example
using OptimizationIpopt

# Declare the ideal solid solution (symbols match symbol(sp) in the database)
ss_ideal    = SolidSolution(["Cal", "Mgs"], IdealSolidSolutionMixingModel())
model_ideal = SolidSolutionActivityModel(DiluteSolutionModel(), [ss_ideal])

opt = IpoptOptimizer(
    acceptable_tol      = 1e-10,
    dual_inf_tol        = 1e-10,
    constr_viol_tol     = 1e-10,
    warm_start_init_point = "no",
)
solver = EquilibriumSolver(cs, model_ideal, opt; variable_space = Val(:log))

# Initial state: 1 L of pure water + 10 mmol of equimolar Ca-Mg carbonate solid
state0 = ChemicalState(cs; T = 298.15u"K", P = 1u"bar")
set_quantity!(state0, "H2O@",  1.0u"kg")
set_quantity!(state0, "Cal",   5e-3u"mol")
set_quantity!(state0, "Mgs",   5e-3u"mol")
set_quantity!(state0, "H+",    1e-7u"mol")
set_quantity!(state0, "OH-",   1e-7u"mol")

state_eq = solve(solver, state0)
nothing #hide
```

```@example ss_example
# Inspect equilibrium amounts [mol]
@printf("%-8s | %12s\n", "Species", "n_eq (mol)")
@printf("%-8s | %12s\n", "--------", "----------")
for sp_i in cs.species
    @printf("%-8s | %12.3e\n", symbol(sp_i), ustrip(us"mol", moles(state_eq, sp_i)))
end
```

### Comparison: pure minerals vs. ideal solid solution

```@example ss_example
# Solver with pure-solid reference (DiluteSolutionModel — no SS)
solver_pure = EquilibriumSolver(cs, DiluteSolutionModel(), opt;
                                variable_space = Val(:log))
state_pure  = solve(solver_pure, state0)
```

```@example ss_example
# Mole fractions of each end-member in the ideal SS at equilibrium
n_cal_eq = ustrip(us"mol", moles(state_eq, "Cal"))
n_mag_eq = ustrip(us"mol", moles(state_eq, "Mgs"))
n_SS     = n_cal_eq + n_mag_eq
x_cal_eq = n_cal_eq / n_SS
x_mag_eq = n_mag_eq / n_SS

@printf("\n%-26s | %s\n", "Quantity", "Value")
@printf("%-26s | %s\n", repeat("-", 26), repeat("-", 12))
@printf("x_Cal at eq.   (SS ideal) | %.4f\n", x_cal_eq)
@printf("x_Mgs at eq.   (SS ideal) | %.4f\n", x_mag_eq)
@printf("log Keff Cal   (SS ideal) | %.2f\n", logKsp_cal - log10(x_cal_eq))
@printf("log Keff Mgs   (SS ideal) | %.2f\n", logKsp_mag - log10(x_mag_eq))

# Dissolved Ca and Mg in the SS vs pure-solid reference
Ca_eq   = ustrip(us"mol", moles(state_eq,   "Ca+2"))
Ca_pure = ustrip(us"mol", moles(state_pure, "Ca+2"))
Mg_eq   = ustrip(us"mol", moles(state_eq,   "Mg+2"))
Mg_pure = ustrip(us"mol", moles(state_pure, "Mg+2"))

@printf("\nn(Ca²⁺) SS ideal  = %.3e mol  (pure Cal = %.3e mol)\n", Ca_eq,   Ca_pure)
@printf("n(Mg²⁺) SS ideal  = %.3e mol  (pure Mgs = %.3e mol)\n", Mg_eq,   Mg_pure)
nothing #hide
```

---

### Non-ideal mixing: the symmetric Margules model

Capobianco & Navrotsky (1987) measured a symmetric interaction parameter
$W \approx 6700$ J/mol for the calcite–magnesite system at 25 °C.

```@example ss_example
W = 6700.0  # J/mol  (Capobianco & Navrotsky 1987)

ss_regular    = SolidSolution(["Cal", "Mgs"], RegularSolidSolutionMixingModel([0.0 W; W 0.0]))
model_regular = SolidSolutionActivityModel(DiluteSolutionModel(), [ss_regular])
solver_reg    = EquilibriumSolver(cs, model_regular, opt; variable_space = Val(:log))
state_reg     = solve(solver_reg, state0)
```

```@example ss_example
n_cal_reg = ustrip(us"mol", moles(state_reg, "Cal"))
n_mag_reg = ustrip(us"mol", moles(state_reg, "Mgs"))
x_cal_reg = n_cal_reg / (n_cal_reg + n_mag_reg)

RT_val  = R * 298.15
lnγ_cal = W * (1 - x_cal_reg)^2 / RT_val
lnγ_mag = W * x_cal_reg^2       / RT_val

@printf("x_Cal     (regular SS) = %.4f\n", x_cal_reg)
@printf("ln γ_Cal               = %+.4f  (>0 → positive deviation)\n", lnγ_cal)
@printf("ln γ_Mgs               = %+.4f\n", lnγ_mag)
nothing #hide
```

!!! note "Positive deviation from ideality"
    The Ca–Mg carbonate system exhibits positive excess Gibbs energy ($W > 0$),
    meaning each end-member is more active than predicted by ideal mixing.
    In practice, this promotes **phase separation** at low temperatures
    (the well-known calcite–magnesite miscibility gap).

---

## API reference

| Type / Function | Description |
|:---|:---|
| [`AbstractSolidSolutionMixingModel`](@ref) | Abstract base for mixing models |
| [`IdealSolidSolutionMixingModel`](@ref) | Ideal mixing: `ln γᵢ = 0` |
| [`RegularSolidSolutionMixingModel`](@ref) | Symmetric Margules: `W` matrix [J/mol] |
| [`lnγ_ss`](@ref) | Compute log activity coefficients |
| [`SolidSolution`](@ref) | Groups end-members + mixing model |
| [`SolidSolutionActivityModel`](@ref) | Decorator wrapping any base model |