# Quickstart

This quickstart shows a few common, minimal examples to get you productive with ChemistryLab in the Julia REPL. It demonstrates parsing formulas, creating species, building reactions and generating a stoichiometric matrix.

## Prerequisites

- Julia 1.6+ (recommended: latest 1.x release)
- Activate the project environment at the repository root (optional for development):

```powershell
# from the project root (PowerShell)
julia --project=.
```

## Install and load

Install the package (or activate the local development copy):

```julia
# In the Julia REPL, enter pkg mode with `]` then:
pkg> add ChemistryLab
```

Then use the package:

```julia
using ChemistryLab
```

## Basic examples

1) Parse a formula and inspect composition

```julia
f = Formula("H2O")
composition(f) # OrderedDict(:H => 2, :O => 1)
charge(f)      # 0
```

2) Create species (with name, symbol and properties)

```julia
water = Species("H2O", name="Water", symbol="H2O", aggregate_state=AS_AQUEOUS)
name(water)    # "Water"
water.M         # molar mass (automatically calculated and stored in properties)
```

3) Build and inspect a Reaction

```julia
r = Reaction("CaCO3 = Ca2+ + CO3-2")
println(r)           # pretty/colored output in the REPL
reactants(r)         # OrderedDict of reactants
products(r)          # OrderedDict of products
```

4) Build a stoichiometric matrix and derive reactions

```julia
# Example species vector
species = [Species("H2O"), Species("H2"), Species("O2")]
A, indep, dep = stoich_matrix(species; display=false)
# Convert a numeric stoichiometric matrix to equation strings
eqns = stoich_matrix_to_equations(A, indep, dep; display=false)
println(eqns)
```

5) Read ThermoFun / Cemdata JSON or .dat files (data directory expected)

```julia
df_elements, df_substances, df_reactions = read_thermofun("data/cemdata18-thermofun.json")
```

## Notes and next steps

- The `Formula`, `Species`, `Reaction` and `stoich_matrix` APIs are intentionally small and composable — explore the `docs/src/` pages for detailed examples.
- For cement-specific workflows, use `CemSpecies` and the `databases` utilities to convert between oxide- and atom-based representations.
- If you want a runnable script that reproduces a specific example from the documentation, I can add one (or convert an example into a Jupyter / Pluto notebook).

Now try the `quickstart` examples interactively in the REPL and then follow the next pages of the tutorial for deeper coverage.
