# Getting started

This quickstart shows a few common, minimal examples to get you productive with ChemistryLab. It demonstrates creating species, building reactions and generating a stoichiometric matrix.

## Simplified example

Let us start with a minimal example in which we compute the thermodynamic properties of a reaction. As a first illustration, we consider the equilibrium of calcite in water. This equilibrium can be written as:

$CaCO_3 \rightleftharpoons Ca^{2+} + {CO_3}^{2-}$

It is possible to calculate the thermodynamic properties of the reaction, in particular the solubility constant of the reaction ($\ln K$) which is related to the Gibbs free energy of the reaction ($\Delta_r G^°$). This solubility constant is a function of temperature and the calculation is performed at a reference temperature of 298 K and at a pressure of 1 Atm using the following equation:

$\Delta_r G^° = - RT\ln K$

where $\Delta_r G^°$ is deduced from the Gibbs energies of formation ($\Delta_f {G_i}^°$) of the other chemical species involved in the reaction:

$\Delta_r G^° = \sum_i \nu_i \Delta_f {G_i}^°$

------------------------

To do this, we can create a list of chemical species, retrieve the thermodynamic properties of these species from one of the databases integrated into ChemistryLab. We can then deduce the chemical species likely to appear in the reaction and calculate the associated stoichiometric matrix.

In this example, the database is [cemdata](https://www.empa.ch/web/s308/thermodynamic-data). The `.json` file is included in ChemistryLab but is a copy of a file which can be found in [ThermoHub](https://github.com/thermohub).

```julia
using ChemistryLab

filebasename = "cemdata18-thermofun.json"
df_elements, df_substances, df_reactions = read_thermofun_database("../../data/" * filebasename)
```

The chemical species likely to appear during calcite equilibrium in water are obtained in the following way:

```@example from_scratch
using ChemistryLab #hide

filebasename = "cemdata18-thermofun.json" #hide
df_elements, df_substances, df_reactions = read_thermofun_database("../../data/" * filebasename) #hide
df_calcite = get_compatible_species(df_substances, split("Cal H2O@ CO2");
                        aggregate_states=[AS_AQUEOUS], exclude_species=split("H2@ O2@ CH4@"), union=true)
dict_species_calcite = build_species_from_database(df_calcite)
```

During species creation, ChemistryLab calculates the molar mass of the species. It also constructs thermodynamic functions (heat capacity, entropy, enthalpy, and Gibbs free energy of formation) as a function of temperature. The evolution of thermodynamic properties as a function of temperature, such as heat capacity, can thus be easily plotted.

```@example from_scratch
dict_species_calcite["Cal"]
```

```@example from_scratch
using Plots

p1 = plot(xlabel="Temperature [K]", ylabel="Cp⁰ [K]", title="Heat capacity of calcite \nas a function of temperature")
plot!(p1, θ -> dict_species_calcite["Cal"].Cp⁰(T = 273.15+θ), 0:0.1:100, label="Cp⁰")
```

Obtaining stoichiometric matrices requires the choice of a species-independent basis.

```@example from_scratch
primaries = [dict_species_calcite[s] for s in split("H2O@ H+ CO3-2 Ca+2")]
SM = StoichMatrix(values(dict_species_calcite), primaries)
```

These stoichiometric matrices thus allow us to write the chemical reactions at work.

```@example from_scratch
list_reactions = reactions(SM)
dict_reactions_calcite = Dict(r.symbol => r for r in list_reactions)
```

Again, when constructing the reactions, the thermodynamic properties of the reactions as a function of temperature are deduced. It is thus possible to see, for example, the expression for the solubility product of calcite for the reaction under study and to plot its evolution.

```@example from_scratch
dict_reactions_calcite["Cal"].logK⁰
```

```@example from_scratch
p2 = plot(xlabel="Temperature [K]", ylabel="pKs", title="Solubility product (pKs) of calcite \nas a function of temperature")
plot!(p2, θ -> dict_reactions_calcite["Cal"].logK⁰(T = 273.15+θ), 0:0.1:100, label="pKs")
```

## Notes and next steps

- The `Formula`, `Species`, `Reaction` and `StoichMatrix` APIs are intentionally small and composable — explore the `docs/src/` pages for detailed examples.
- For cement-specific workflows, use `CemSpecies` and the `databases` utilities to convert between oxide- and atom-based representations.

Now try the `quickstart` examples interactively in the REPL and then follow the next pages of the tutorial for deeper coverage.