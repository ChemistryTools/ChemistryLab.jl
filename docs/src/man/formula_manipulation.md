# Chemical Formula Manipulation

ChemistryLab allows you to create and manipulate chemical formulas. It provides the `Formula` type (holds a formula string, multiple string representations, an atom composition map and a charge) and `AtomGroup` helper values, plus utilities to convert between notations (Phreeqc ↔ Unicode), to format/colour formulas, to perform arithmetic on formulas, and to validate atom symbols.

```julia
struct Formula{T<:Number}
    expr::String
    phreeqc::String
    unicode::String
    colored::String
    composition::OrderedDict{Symbol,T}
    charge::Int8
end
```

## Formula construction

 Formulas can be constructed:
- by parsing a string containing eventually fractional or decimal coefficients
```@example 1
using ChemistryLab #hide
f = Formula("C3AFS5//8H4.32")
```

- from a dictionary
```@example 1
using ChemistryLab #hide
fCaCO3 = Formula(Dict(:Ca => 1, :C => 1, :O => 3))
```

- from atom groups
```@example
using ChemistryLab #hide
fCO2 = AtomGroup(:C) + AtomGroup(2,:O)
```

Charges can also be included during the creation in two different ways:
```@example 1
fHSO₄⁻ = AtomGroup(:H)+AtomGroup(:S)+AtomGroup(4,:O)+AtomGroup(:e)
```

Or:
```@example 1
fNa⁺ = AtomGroup(:Na)+AtomGroup(:Zz)
```

## Type of Formula

The type of the `Formula` `struct` being associated with the most complex type of the set of coefficients.

```@example 1
typeof(Formula("H2O"))
```

```@example 1
typeof(Formula("C3AFS5//8H4.32"))
```

## Change of type

Coefficient types can be converted *a posteriori*.

```@example 1
convert(Float64, Formula("H2O"))
```



---


