# Equilibrium

```@index
Pages = ["equilibrium.md"]
```

## Activity models

```@docs
AbstractActivityModel
DiluteSolutionModel
HKFActivityModel
DaviesActivityModel
activity_model
build_potentials
REJ_HKF
REJ_CHARGE_DEFAULT
hkf_debye_huckel_params
```

## Solid solutions

```@docs
AbstractSolidSolutionModel
IdealSolidSolutionModel
RedlichKisterModel
AbstractSolidSolutionPhase
SolidSolutionPhase
end_members
model
```

## Problem and solver

```@docs
EquilibriumProblem
EquilibriumSolver
equilibrate
```
