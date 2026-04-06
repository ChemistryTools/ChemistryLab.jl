# Chemical systems and states

```@index
Pages = ["chemical_systems.md"]
```

## ChemicalSystem

```@docs
ChemicalSystem
aqueous
crystal
gas
solutes
solvent
components(cs::ChemicalSystem)
gasfluid
get_reaction
```

## ChemicalState

```@docs
ChemicalState
temperature
pressure
set_temperature!
set_pressure!
moles
set_quantity!
rescale!
mass
volume
pH
pOH
porosity
saturation
```
