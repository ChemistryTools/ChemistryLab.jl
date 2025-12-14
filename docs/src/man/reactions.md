# Chemical reactions

In ChemistryLab it is possible to build chemical reactions and manipulate them. A reaction is constructed as a structure, "a composite data type that allows you to store multiple values in a single object". The `struct` is organized as follows:

```julia
struct Reaction{SR<:AbstractSpecies, TR<:Number, SP<:AbstractSpecies, TP<:Number}
    equation::String
    colored::String
    reactants::OrderedDict{SR, TR}
    products::OrderedDict{SP, TP}
    equal_sign::Char
    properties::OrderedDict{Symbol,PropertyType}
end
```

## Parsing reactions

`Reaction` is a composite type `struct` which can be build from:

- a string containing [species](./databases.md#species)

```julia
equation = "13H⁺ + NO₃⁻ + CO₃²⁻ + 10e⁻ = 6H₂O@ + HCN@"
reac, prod, equal_sign = parse_equation(equation)
```

- a string containing [cement species](./databases.md#species)

```julia
eqC3S = "C₃S + 5.3H = 1.3CH + C₁.₇SH₄"
rC3S = CemReaction(eqC3S)
```

- an operation on species

```@example
using ChemistryLab
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
CSH = CemSpecies("C1.7SH4")
r = C3S + 5.3H ↔ 1.3CH + CSH
typeof(r)
```

- a balance calculation

```@example
using ChemistryLab
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
CSH = CemSpecies("C1.7SH4")
r = Reaction([C3S, H, CH, CSH]; equal_sign='→')
pprint(r)
```

- a balance calculation with symbolic numbers

```julia
using ChemistryLab
using SymPy
a, b, g = symbols("a b g", real=true)
CSH = CemSpecies(Dict(:C => a, :S => one(Sym), :H => g))
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
r = Reaction([CSH, C3S, H, CH]; equal_sign='→')
```

```julia
using ChemistryLab
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
CSH = CemSpecies("C1.7SH4")
r = map(simplify, Reaction([C3S, H], [CH, CSH]; equal_sign='→'))
SM = StoichMatrix([C3S], [CSH, H, CH])
```

!!! note "Collection of arrow symbols"
    In ChemistryLab, there are collections of arrow symbols used in chemical reaction notation, such as:
    - `>`, `→`, `↣`, `↦`, `⇾`, `⟶`, `⟼`, `⥟`, `⥟`, `⇀`, `⇁`, `⇒`, `⟾` for reaction directionality from reactants to products;
    - `<`, `←`, `↢`, `↤`, `⇽`, `⟵`, `⟻`, `⥚`, `⥞`, `↼`, `↽`, `⇐`, `⟽` for reaction directionality from products to reactants;
    - `↔`, `⟷`, `⇄`, `⇆`, `⇌`, `⇋`, `⇔`, `⟺` for reversible reactions and equilibrium states;
    - `=`, `≔`, `⩴`, `≕` to separate reactants from products in balanced equations.