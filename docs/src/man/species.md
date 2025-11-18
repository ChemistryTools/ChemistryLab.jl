
# Species

`Species` is a composite type (introduced by the keyword `struct`) and is defined by a human-readable name, a chemical symbol/notation, an underlying `Formula` object holding composition, charge, and string representations, a physical state `aggregate_state`, a species `class`, as well as an extensible map of custom properties (molar mass, thermodynamic data, etc.)


```julia
struct Species{T<:Number} <: AbstractSpecies
    name::String
    symbol::String
    formula::Formula{T}
    aggregate_state::AggregateState
    class::Class
    properties::OrderedDict{Symbol,PropertyType}
end
```

!!! info "Advanced description"
    - `aggregate_state` denotes the state of the species (solid, liquid, gas) for which the possible keywords are `AS_AQUEOUS`, `AS_CRYSTAL`, `AS_GAS` and `AS_UNDEF`
    - `class` defines the role played by the species in the solution. The possible keywords are `SC_AQSOLVENT`, `SC_AQSOLUTE`, `SC_COMPONENT`, `SC_GASFLUID` and `SC_UNDEF`
    - `properties` refers to the set of properties intrinsic to the species. These properties are detailed below. 

## Species construction

`Species` can be created from:
- a `Formula`

```@example
using ChemistryLab #hide
fH2O = Formula("H2O")
H2O = Species(fH2O)
```

- a string
```@example
using ChemistryLab #hide
HSO4⁻ = Species("HSO₄⁻")
```

- a dictionary
```@example
using ChemistryLab #hide
CO2 = Species(Dict(:C => 1, :O => 2))
```

!!! note "Adding charge"
    To add a charge when creating species with a dictionary, you must add, after the dictionary, the value of the charge (charge is considered an argument of the composite type).

```@example
using ChemistryLab #hide
CO2 = Species(Dict(:Si => 1, :O => 3), -2)
```

Keyword arguments such as `name`, `symbol`, `aggregate_state`, `class` can be added during construction.

```@example H2O
using ChemistryLab #hide
fH₂O = Formula("H2O")
H₂O = Species(fH₂O; name="Water", symbol="H₂O@", aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)
```

And `symbol` accept unicode characters.
```@example CO2
using ChemistryLab #hide
CO₂ = Species(Dict(:C=>1, :O=>2); name="Carbon dioxide", symbol="CO₂⤴", aggregate_state=AS_GAS, class=SC_GASFLUID)
```

!!! note "Comparison between species"
    Comparison between species (or cemspecies) are done by comparing atoms, aggregate_state and class. In the example below, vapour is not equal to H₂O since *aggregate_state* and *class* are different despite atoms are identical.
    ```julia
    vapour = Species("H2O"; name="Vapour", symbol="H₂O⤴", aggregate_state=AS_GAS, class=SC_GASFLUID)
    vapour == H₂O
    ```


!!! tip "Remark"
    You will also have noticed that a calculation of the molar mass of the species is systematically carried out.

---

## Species properties

The molar mass is systematically calculated and integrated into the species properties. The heat capacity function is also integrated as a predefined function of the `Species` structure. This function is expressed as follows:

${C_p}^° = a_0 + a_1 T + a_2 T^{−2} + a_3 T^{−0.5}$

Other species properties are open and left to the discretion of users. We can of course imagine that these properties could contain thermodynamic properties such as the Gibbs energy of formation or even the entropy variation, these properties themselves being temperature dependent. These properties must nevertheless respect one of the following types: `Number`, `AbstractVector{<:Number}`, `Function`, `AbstractString`.

Imagine, for example, that we wanted to construct the $CO_2$ molecule with some of its thermodynamic properties. The Gibbs energy of formation of this species is equal to $-394.39\; KJ/mol$. This property, intrinsic to the species, can be added simply as follows:

```@example CO2
using DynamicQuantities, ModelingToolkit
CO₂.ΔfG⁰ = -394.39u"kJ/mol"
```

Heat capacity, on the other hand, is introduced in the following way: 

```@example CO2
coeffs = [:a₀ => 44.22u"J/K/mol", :a₁ => 0.0088u"J/mol/K^2", :a₂ => -861.904e6u"J*K/mol", :a₃ => 0.0u"J/mol/√K"]
CO₂.Cp = ThermoFunction(:Cp, coeffs; ref=[:T=>298.15u"K", :P=>1u"bar"])
```

!!! note "Heat capacity value"
    By default, the call of $C_p$ function return the value of the function for a temperature equal to $T_{ref}$.
    ```julia
    CO₂.Cp()
    ```