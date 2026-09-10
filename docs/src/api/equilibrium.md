# Equilibrium

```@index
Pages = ["equilibrium.md"]
```

## Activity models

```@autodocs
Modules = [ChemistryLab]
Pages = ["equilibrium/activities.jl"]
```

## Aqueous properties

Molalities, ionic strength, activity coefficients and the activity-convention
pH, read back off a solved state. See
[Reading the aqueous properties back](@ref sec-aqueous-properties) for the two
traps these exist to avoid.

```@autodocs
Modules = [ChemistryLab]
Pages = ["equilibrium/aqueous_properties.jl"]
```

## Solid solutions

```@autodocs
Modules = [ChemistryLab]
Pages = ["chemical_structs/solid_solutions.jl"]
```

## Certified equilibrium

Newton on the KKT system in element-potential space, and the certificate that
proves a composition optimal. See
[Proving that an answer is the answer](@ref sec-theory-certificate) for the derivation, and
[`equilibrate_certified`](@ref) for the route `equilibrate` takes by default.

```@autodocs
Modules = [ChemistryLab]
Pages = ["equilibrium/dual_solver.jl", "equilibrium/certified.jl"]
```

## Constraints

What is held fixed while the Gibbs energy is minimized. See
[Constraints other than fixed T and P](@ref sec-equilibrium-constraints).

```@autodocs
Modules = [ChemistryLab]
Pages = ["equilibrium/constraints.jl"]
```

## Water retention

The relation between how much water a pore space still holds and how tightly it
holds it — the constitutive input [`CapillaryWater`](@ref) needs, and the source
of the humidity [`PoreHumidity`](@ref) hands to a rate law. It is measured, not
assumed, so nothing here carries a default value.

```@autodocs
Modules = [ChemistryLab]
Pages = ["equilibrium/retention.jl"]
```

## Problem and solver

```@autodocs
Modules = [ChemistryLab]
Pages = ["equilibrium/equilibrium_problems.jl", "equilibrium/equilibrium_solver.jl"]
```
