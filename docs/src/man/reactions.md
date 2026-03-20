# Chemical Reactions

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
C3S = CemSpecies("C3S", symbol="C₃S", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
H = CemSpecies("H", symbol="H₂O@", aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)
CH = CemSpecies("CH", symbol="C₃S", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
CSH = CemSpecies("C1.7SH4", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
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

```@example
using ChemistryLab
using ModelingToolkit
@variables a b g
CSH = CemSpecies(Dict(:C => a, :S => one(Num), :H => g))
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
r = Reaction([CSH, C3S, H, CH]; equal_sign='→')
```

```@example
using ChemistryLab
using ModelingToolkit
using PrettyTables
@variables a b g
CSH = CemSpecies(Dict(:C => a, :S => one(Num), :H => g))
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
r = map(simplify, Reaction([C3S, H], [CH, CSH]; equal_sign='→'))
SM = StoichMatrix([C3S], [CSH, H, CH])
pprint(SM)
```

!!! note "Collection of arrow symbols"
    In ChemistryLab, there are collections of arrow symbols used in chemical reaction notation, such as:
    - `>`, `→`, `↣`, `↦`, `⇾`, `⟶`, `⟼`, `⥟`, `⥟`, `⇀`, `⇁`, `⇒`, `⟾` for reaction directionality from reactants to products;
    - `<`, `←`, `↢`, `↤`, `⇽`, `⟵`, `⟻`, `⥚`, `⥞`, `↼`, `↽`, `⇐`, `⟽` for reaction directionality from products to reactants;
    - `↔`, `⟷`, `⇄`, `⇆`, `⇌`, `⇋`, `⇔`, `⟺` for reversible reactions and equilibrium states;
    - `=`, `≔`, `⩴`, `≕` to separate reactants from products in balanced equations.

## Thermodynamic properties of reactions

When the species involved in a reaction carry thermodynamic data (loaded from a database), the reaction automatically exposes temperature-dependent thermodynamic functions. These are computed lazily on first access and stored in the reaction's `properties` dict:

| Property | Description |
|----------|-------------|
| `r.ΔᵣCp⁰` | Heat capacity of reaction (J mol⁻¹ K⁻¹) |
| `r.ΔᵣH⁰`  | Enthalpy of reaction (J mol⁻¹) |
| `r.ΔᵣS⁰`  | Entropy of reaction (J mol⁻¹ K⁻¹) |
| `r.ΔᵣG⁰`  | Gibbs free energy of reaction (J mol⁻¹) |
| `r.logK⁰` | Decimal logarithm of the equilibrium constant |

Each property is a `SymbolicFunc` callable with a keyword argument `T` (temperature in K):

```julia
using ChemistryLab

# Load reactions from a database-built stoichiometric matrix
all_species = build_species("../../../data/cemdata18-merged.json")
species = speciation(all_species, split("Cal H2O@");
              aggregate_state=[AS_AQUEOUS], exclude_species=split("H2@ O2@ CH4@"))
dict_species = Dict(symbol(s) => s for s in species)
candidate_primaries = [dict_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_species, s)]
cs = ChemicalSystem(species, candidate_primaries)
list_reactions = reactions(cs.SM)
dict_reactions = Dict(r.symbol => r for r in list_reactions)

r_cal = dict_reactions["Cal"]   # calcite dissolution reaction

# Evaluate at 25 °C
r_cal.logK⁰(T = 298.15)          # log₁₀ K at 25 °C
r_cal.ΔᵣG⁰(T = 298.15)           # ΔᵣG° at 25 °C  (J/mol)
```

Properties can also be set manually on a reaction:

```julia
r = Reaction("H2 + O2 = H2O")
r[:ΔᵣH⁰] = -241800.0   # J/mol
```

## Reactions from a stoichiometric matrix

The most common source of reactions in a database workflow is `reactions(SM)`, which derives all independent reactions from a `StoichMatrix`:

```julia
list_reactions = reactions(cs.SM)
dict_reactions = Dict(r.symbol => r for r in list_reactions)
```

Each reaction symbol matches the dependent species it describes. Temperature-dependent thermodynamic functions are computed automatically when the constituent species carry thermodynamic data.

---