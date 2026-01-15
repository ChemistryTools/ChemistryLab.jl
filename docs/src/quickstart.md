# Quickstart

This quickstart shows a few common, minimal examples to get you productive with ChemistryLab in the Julia REPL. It demonstrates parsing formulas, creating species, building reactions and generating a stoichiometric matrix.

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

- Parse a formula and inspect composition

```julia
f = Formula("H2O")
composition(f) # OrderedDict(:H => 2, :O => 1)
charge(f)      # 0
```

- Create species (with name, symbol and properties)

```julia
water = Species("H2O", name="Water", symbol="H2O", aggregate_state=AS_AQUEOUS)
name(water)    # "Water"
water.M         # molar mass (automatically calculated and stored in properties)
```

- Build and inspect a Reaction

```julia
r = Reaction("CaCO3 = Ca+2 + CO32-")
println(r)           # pretty/colored output in the REPL
reactants(r)         # OrderedDict of reactants
products(r)          # OrderedDict of products
```

- Build a stoichiometric matrix and derive reactions

```julia
# Example species vector
species = Species.(split("H₂O H⁺ OH⁻ CO₂ HCO₃⁻ CO₃²⁻"))
SM = StoichMatrix(species)
# Convert a numeric stoichiometric matrix to independent reactions
list_reactions = reactions(SM)
pprint(list_reactions)
```

- Read ThermoFun JSON files (data directory expected)

```julia
df_elements, df_substances, df_reactions, dict_species, dict_reactions = read_thermofun("data/cemdata18-thermofun")
```

## Notes and next steps

- The `Formula`, `Species`, `Reaction` and `StoichMatrix` APIs are intentionally small and composable — explore the `docs/src/` pages for detailed examples.
- For cement-specific workflows, use `CemSpecies` and the `databases` utilities to convert between oxide- and atom-based representations.

Now try the `quickstart` examples interactively in the REPL and then follow the next pages of the tutorial for deeper coverage.
