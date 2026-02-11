
# [Species](@id sec-species)

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

$a_0 + a_1 * T + a_2 * T^{-2} + a_3 * T^{-0.5} + a_4 * T^2 + a_5 * T^3 + a_6 * T^4 + a_7 * T^{-3} + a_8 * T^{-1} + a_9 * T^{0.5} + a_{10} * log(T)$

Other species properties are open and left to the discretion of users. We can of course imagine that these properties could contain thermodynamic properties such as the Gibbs energy of formation or even the entropy variation, these properties themselves being temperature dependent. These properties must nevertheless respect one of the following types: `Number`, `AbstractVector{<:Number}`, `Function`, `AbstractString`.

Imagine, for example, that we wanted to construct the $CO_2$ molecule in a gaseous state with some of its thermodynamic properties. The thermodynamic properties of the molecule, which can be found for example on the website [thermoddem](https://thermoddem.brgm.fr/), are as follows at 298 K and 1 atm:

- heat capacity ${C_p}^°$: 37.14 $J.mol^{-1}.K^{-1}$
- molar volume $V^°$: 25.3 $cm^3.mol^{-1}$
- enthalpy of formation $\Delta_a {H^°}$: -393510 $J.mol^{-1}$
- entropy $S^°$: 213.785 $J.mol^{-1}.K^{-1}$ 
- Gibbs free energy of formation $\Delta_a {G^°}: -394373 $J.mol^{-1}$

Furthermore, as explained above, heat capacity is a function of temperature. The parameters $a_0$, $a_1$, $a_2$, and $a_3$ can also be found on the same website. For CO2, the values ​​are as follows: $a_0 = 33.98$, $a_1 = 23.88e-3$, $a_2 = 0$ et $a_3 = 0$. 

```julia
using DynamicQuantities, ModelingToolkit
th_prop_0_CO2 = Dict(:Cp⁰ => 37.14, :ΔₐH⁰ => -393510, :S⁰ => 213.785, :ΔₐG⁰ => -394373, :V⁰ => 25.3)
coeffs = Dict(:a₀ => 33.98, :a₁ => 23.88e-3, :a₂ => 0.0, :a₃ => 0.0)
```

!!! note "Heat capacity function"
    Although the function describing heat capacity has many parameters, it is of course possible to use only some of them. Here, only the parameters $a_0$ and $a_1$ are non-zero. The expression therefore becomes: $a_0 + a_1 * T$.


!!! danger "Reference temperature"
    It is important to define the reference temperature at which the thermodynamic properties are measured.
    ```julia
    T_ref = Dict(:T => 298.15)
    ```

#### Heat capacity, enthalpy and free energy as a function of temperature

Reference thermodynamical properties and temperature being defined, a simple call to `build_thermo_functions` allows the thermodynamic functions, such as heat capacity, entropy, enthalpy, and free enthalpy to be built as a function of temperature.

```@example CO2
using DynamicQuantities, ModelingToolkit #hide
th_prop_0_CO2 = Dict(:Cp⁰ => 37.14, :ΔₐH⁰ => -393510, :S⁰ => 213.785, :ΔₐG⁰ => -394373, :V⁰ => 25.3) #hide
coeffs = Dict(:a₀ => 33.98, :a₁ => 23.88e-3, :a₂ => 0.0, :a₃ => 0.0) #hide
T_ref = Dict(:T => 298.15) #hide
params_Cp_CO2 = merge(th_prop_0_CO2, coeffs, T_ref)
dtf_CO2 = build_thermo_functions(:cp_ft_equation, params_Cp_CO2)
```
