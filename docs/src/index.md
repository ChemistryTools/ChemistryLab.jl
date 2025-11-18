```@meta
CurrentModule = ChemistryLab
```

# ChemistryLab

Documentation for [ChemistryLab](https://github.com/jfbarthelemy/ChemistryLab.jl).

ChemistryLab is a Julia package for creating and manipulating chemical formulas and species, building stoichiometric matrices and reactions, and importing thermodynamic data (ThermoFun / Cemdata) with callable thermodynamic functions.

The library aims to model the equilibrium of geochemical systems by minimizing the Gibbs energy, as well as offering flexible and extensible modeling capabilities for chemical equilibrium and kinetics calculations.

There are many applications such as mineral dissolution and precipitation, aqueous speciation, but also cement hydration.

## Install

```julia-repl
pkg> add https://github.com/jfbarthelemy/ChemistryLab.jl
```



## Quickstart

Install ChemistryLab in your chosen environment by entering pkg mode by pressing `]` and then:

```julia
pkg> add ChemistryLab
```

In order to use ChemistryLab, it is then necessary to load the ChemistryLab.jl package:

```julia
julia> using ChemistryLab
```
