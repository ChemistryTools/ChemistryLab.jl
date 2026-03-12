# [Chemical Equilibrium](@id sec-equilibrium)

ChemistryLab computes thermodynamic equilibrium by minimising the Gibbs free energy of the system subject to element-conservation constraints. The workflow always follows the same four steps:

1. Build a [`ChemicalSystem`](@ref) (species + stoichiometric matrix).
2. Create an initial [`ChemicalState`](@ref) (temperature, pressure, initial amounts).
3. Call [`equilibrate`](@ref) (or use [`EquilibriumSolver`](@ref) explicitly).
4. Inspect the resulting [`ChemicalState`](@ref).

---

## Minimal workflow

The convenience function [`equilibrate`](@ref) handles everything with sensible defaults.
The example below computes the equilibrium state of calcite (CaCO₃) dissolving in mildly acidic water — a standard geochemical benchmark.

```@example eq_setup
using ChemistryLab
using DynamicQuantities

substances = build_species("../../../data/slop98-inorganic-thermofun.json")

# Select the carbonate-system species, calcite and its dissolution product Ca²⁺
dict = Dict(symbol(s) => s for s in substances)
species = [dict[sym] for sym in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")]

cs = ChemicalSystem(species, ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
```

```@example eq_setup
state = ChemicalState(cs)

# 1 mmol calcite dissolved in 1 L of acidic water (initial pH ≈ 4)
set_quantity!(state, "Cal",  1e-3u"mol")
set_quantity!(state, "H2O@", 1.0u"kg")

V = volume(state)
set_quantity!(state, "H+",  1e-4u"mol/L" * V.liquid)   # pH = 4
set_quantity!(state, "OH-", 1e-10u"mol/L" * V.liquid)  # charge seed

state_eq = equilibrate(state)
```

!!! tip "Quick shortcut"
    Calling `equilibrate(state)` with no extra arguments uses sensible defaults and is usually sufficient for aqueous geochemical problems.

---

## Inspecting the equilibrium state

The returned [`ChemicalState`](@ref) carries all derived thermodynamic quantities:

```@example eq_setup
println("pH      = ", pH(state_eq))
println("pOH     = ", pOH(state_eq))
println("porosity   = ", porosity(state_eq))
println("saturation = ", saturation(state_eq))
```

Phase volumes and mole amounts are accessible via named tuples:

```@example eq_setup
v = volume(state_eq)
println("V liquid = ", v.liquid)
println("V solid  = ", v.solid)
println("V total  = ", v.total)

m = moles(state_eq)
println("n liquid = ", m.liquid)
println("n solid  = ", m.solid)
```

Individual species amounts (in mol):

```@example eq_setup
cs_eq = state_eq.system
for (i, sp) in enumerate(cs_eq.species)
    n_i = state_eq.n[i]
    println(rpad(symbol(sp), 20), ustrip(n_i), " mol")
end
```

---

## Controlling the solver

### Variable space: `:linear` vs `:log`

[`equilibrate`](@ref) accepts a `variable_space` keyword that selects the optimisation variable space:

| `variable_space`        | Variables | Recommended when |
|------------------|-----------|-----------------|
| `Val(:linear)`   | mole amounts `nᵢ ≥ 0` | most systems, default |
| `Val(:log)`      | `log nᵢ` | systems spanning many orders of magnitude |

```julia
state_eq_log = equilibrate(state; variable_space=Val(:log))
```

!!! warning "Convergence"
    Solving a system of equations in chemistry can be a difficult undertaking. The orders of magnitude can vary greatly, and convergence is not guaranteed.

---

### Tolerances

Tighter tolerances are passed directly as keyword arguments and forwarded to the underlying Ipopt solver:

```julia
state_eq_tight = equilibrate(state; abstol=1e-12, reltol=1e-12)
```

---

## Using `EquilibriumSolver` explicitly

For batch calculations where many different initial states share the same system and activity model, construct an [`EquilibriumSolver`](@ref) once and reuse it:

```julia
using OptimizationIpopt

opt = IpoptOptimizer(
    acceptable_tol        = 1e-12,
    dual_inf_tol          = 1e-12,
    acceptable_iter       = 1000,
    constr_viol_tol       = 1e-12,
    warm_start_init_point = "no",
)

solver = EquilibriumSolver(
    cs,
    DiluteSolutionModel(),
    opt;
    variable_space = Val(:linear),
    abstol  = 1e-10,
    reltol  = 1e-10,
)
```

Once built, `solver` is called with any compatible `ChemicalState`:

```@example eq_setup
using OptimizationIpopt #hide

opt = IpoptOptimizer( #hide
    acceptable_tol        = 1e-12, #hide
    dual_inf_tol          = 1e-12, #hide
    acceptable_iter       = 1000, #hide
    constr_viol_tol       = 1e-12, #hide
    warm_start_init_point = "no", #hide
) #hide

solver = EquilibriumSolver( #hide
    cs, #hide
    DiluteSolutionModel(), #hide
    opt; #hide
    variable_space = Val(:linear), #hide
    abstol  = 1e-10, #hide
    reltol  = 1e-10, #hide
) #hide
state_eq2 = solve(solver, state)
```

!!! note "Performance"
    The potential function `μ(n, p)` is compiled once during `EquilibriumSolver` construction. Repeated calls to `solve(solver, ...)` with different states reuse it, avoiding redundant compilation overhead.

---

## Temperature dependence (10–30 °C)

Calcite solubility varies with temperature. Using the `solver` built above, we sweep from 10 to 30 °C and track pH, dissolved calcium and remaining solid calcite:

```@example eq_setup
using Plots


temperatures = 10:30   # °C

pH_vals   = Float64[]
nCa_vals  = Float64[]  # mmol
nCal_vals = Float64[]  # mmol

i_Ca  = findfirst(sp -> symbol(sp) == "Ca+2", cs.species)
i_Cal = findfirst(sp -> symbol(sp) == "Cal",  cs.species)

s = ChemicalState(cs)
for θ in temperatures
    set_temperature!(s, (273.15 + θ) * u"K")
    s_eq = solve(solver, s)
    push!(pH_vals,   pH(s_eq))
    push!(nCa_vals,  ustrip(s_eq.n[i_Ca]) * 1e3)
    push!(nCal_vals, ustrip(s_eq.n[i_Cal]) * 1e3)
end
```

The figures can then be drawn.

```@example eq_setup
p1 = plot(collect(temperatures), pH_vals,
    xlabel = "T (°C)", ylabel = "pH", label = "pH",
    marker = :circle, linewidth = 2, title = "pH")
p2 = plot(collect(temperatures), nCa_vals,
    xlabel = "T (°C)", ylabel = "n (mmol)", label = "Ca²⁺",
    marker = :circle, linewidth = 2, title = "Dissolved species")
plot!(p2, collect(temperatures), nCal_vals,
    label = "Cal", marker = :square, linewidth = 2)
plot(p1, p2, layout = (1, 2), size = (700, 350))
```

!!! note "Calcite solubility"
    Calcite is a **retrograde soluble** mineral: its solubility decreases with increasing temperature, so less Ca²⁺ is released and pH rises slightly as temperature increases.

---

## Activity models

All activity models inherit from [`AbstractActivityModel`](@ref). The only built-in model is [`DiluteSolutionModel`](@ref), which implements:

| Phase | Law | Expression |
|-------|-----|-----------|
| Solvent (H₂O) | Raoult | `ln a = ln xₛ` |
| Aqueous solutes | Henry | `ln a = ln(cᵢ / c°)`, `c° = 1 mol/L` |
| Crystals | Pure solid | `ln a = 0` |
| Gas | Ideal mixture | `ln a = ln xᵢ` |

To implement a custom activity model, define a new subtype and extend `activity_model`:

```julia
struct MyModel <: AbstractActivityModel
    # model parameters
end

function ChemistryLab.activity_model(::ChemicalSystem, ::MyModel)
    # return a closure lna(n, p) -> Vector{Float64}
    return (_, _) -> begin
        # compute and return log-activities
    end
end
```

Pass your model to `equilibrate` or `EquilibriumSolver`:

```julia
state_eq = equilibrate(state; model=MyModel(...))
```
