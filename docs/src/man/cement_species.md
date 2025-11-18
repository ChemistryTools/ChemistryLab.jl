

# Cement Species

The manipulation of chemical formulas can also be done in cement notation. `CemSpecies` is a composite type similar to `Species`.


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

## CemSpecies construction

Here are several example of `CemSpecies` construction:

```@setup example_cemspecies
    using ChemistryLab
```

```@example example_cemspecies
C3S = CemSpecies("C3S"; name="Alite", symbol="C₃S", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
C2S = CemSpecies("C₂S"; name="Belite", symbol="C₂S", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
C3A = CemSpecies("C3A"; name="Aluminate", symbol="C₃A", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
C4AF = CemSpecies(Dict(:C => 4, :A => 1, :F => 1); name="Ferrite", symbol="C₄AF")
```

!!! warning "Warning"
    Not every molecule can be used to build a cement species. It is necessary for this molecule to decompose into a combination of the oxides present in the manufacturers' cement sheet (e.g. $CaO$, $SiO_2$, $Fe_2O_3$, $Al_2O_3$) and water. Thus, the following code will return an error.
    ```julia
    CemSpecies(Species("Ca(OH)"))
    ```


---

## Numeric and Symbolic CemSpecies

The previous species were constructed from integer values ​​of the number of chemical elements. However, other numerical value types ​​are possible (as for [formulas](./formula_manipulation.md#formulas)), such as fraction or Real values.

```@example
using ChemistryLab
ox = Dict(:C => 1.666667, :S => 1, :H => 2.1)
jennite = CemSpecies(ox)
```

Symbolic values are also allowed. In this case, you need to use the [`SymPy`](https://github.com/JuliaPy/SymPy.jl) library:

```@example sympy1
using ChemistryLab
using SymPy
â, ĝ = symbols("â ĝ", real = true)
ox = Dict(:C => â, :S => one(Sym), :H => ĝ)
CSH = CemSpecies(ox)
```

The value of variables can be defined *a posteriori*.

```julia
jennite = CemSpecies(map(N, map(subs, cemformula(CSH), â => 1.666667, ĝ => 2.1)))
```

!!! note "Remark"
    Conversion of coefficient types can also be done.
    ```julia
    floatCSH = Species(convert(Float64, formula(numCSH)))
    ```

---

## Conversion of Species to CemSpecies and vice versa

Convert species to cement notation and Unicode. Conversion can be done on simple species:

```@example example_cemspecies
H2O = Species("H₂O")
cemH2O = CemSpecies(H2O)
```

```@example example_cemspecies
C3S = CemSpecies("C3S"; name="Alite", symbol="C₃S", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
spC3S = Species(C3S)
```

Or more complex one:

```@example CSH
using ChemistryLab #hide
CSH = Species("(SiO2)1(CaO)1.666667(H2O)2.1")
jennite = CemSpecies(CSH)
```

---
