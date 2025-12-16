<p>
  <img src="./docs/src/assets/logo.svg" width="100">
</p>

# ChemistryLab

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jfbarthelemy.github.io/ChemistryLab.jl/dev/)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jfbarthelemy.github.io/ChemistryLab.jl/stable/)
[![Build Status](https://github.com/jfbarthelemy/ChemistryLab.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jfbarthelemy/ChemistryLab.jl/actions/workflows/CI.yml?query=branch%3Amain)

[![DOI](https://zenodo.org/badge/1054296488.svg)](https://doi.org/10.5281/zenodo.17756074)

ChemistryLab.jl is a computational chemistry toolkit. Although initially dedicated to low-carbon cementitious materials and aqueous solutions and designed for researchers, engineers, and developers working with cement chemistry, its scope is actually wider. It provides formula handling, species management, stoichiometric matrix construction, and database interoperability (ThermoFun and Cemdata). Main features include chemical formula parsing, Unicode/Phreeqc notation conversion, reaction and equilibrium analysis, and data import/export.

## Features

- **Chemical formula handling**: Create, convert, and display formulas with charge management and Unicode/Phreeqc notation.
- **Chemical species management**: `Species` and `CemSpecies` types to represent solution and solid phase species.
- **Stoichiometric matrices**: Automatic construction of matrices for reaction and equilibrium analysis.
- **Database interoperability**: Import and merge ThermoFun (.json) and Cemdata (.dat) data.
- **Parsing tools**: Convert chemical notations, extract charges, calculate molar mass, and more.

## Examples

- Reaction defined from a string

```julia
julia> using ChemistryLab

julia> equation = "13H⁺ + NO₃⁻ + CO₃²⁻ + 10e⁻ = 6H₂O@ + HCN@"
"13H⁺ + NO₃⁻ + CO₃²⁻ + 10e⁻ = 6H₂O@ + HCN@"

julia> r = Reaction(equation)
13H⁺ + NO₃⁻ + CO₃²⁻ + 10e⁻ = 6H₂O@ + HCN@
 reactants: H⁺ => 13, NO₃⁻ => 1, CO₃²⁻ => 1
  products: H₂O@ => 6, HCN@ => 1
    charge: -10
```

- Stoichiometric matrix construction

```julia
julia> H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻ = Species.(split("H₂O H⁺ OH⁻ CO₂ HCO₃⁻ CO₃²⁻"))
6-element Vector{Species{Int64}}:
 H₂O {H₂O} [H₂O ◆ H2O]
 H⁺ {H⁺} [H⁺ ◆ H+]
 OH⁻ {OH⁻} [OH⁻ ◆ OH-]
 CO₂ {CO₂} [CO₂ ◆ CO2]
 HCO₃⁻ {HCO₃⁻} [HCO₃⁻ ◆ HCO3-]
 CO₃²⁻ {CO₃²⁻} [CO₃²⁻ ◆ CO3-2]

julia> CSM = CanonicalStoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻])
┌────┬─────┬────┬─────┬─────┬───────┬───────┐
│    │ H₂O │ H⁺ │ OH⁻ │ CO₂ │ HCO₃⁻ │ CO₃²⁻ │
├────┼─────┼────┼─────┼─────┼───────┼───────┤
│  C │     │    │     │   1 │     1 │     1 │
│  H │   2 │  1 │   1 │     │     1 │       │
│  O │   1 │    │   1 │   2 │     3 │     3 │
│ Zz │     │  1 │  -1 │     │    -1 │    -2 │
└────┴─────┴────┴─────┴─────┴───────┴───────┘

julia> SM = StoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻]) # selection of a basis of independent primaries
┌─────┬─────┬────┬─────┬─────┬───────┬───────┐
│     │ H₂O │ H⁺ │ OH⁻ │ CO₂ │ HCO₃⁻ │ CO₃²⁻ │
├─────┼─────┼────┼─────┼─────┼───────┼───────┤
│ H₂O │   1 │    │   1 │     │     1 │     1 │
│  H⁺ │     │  1 │  -1 │     │    -1 │    -2 │
│ CO₂ │     │    │     │   1 │     1 │     1 │
└─────┴─────┴────┴─────┴─────┴───────┴───────┘

julia> reactions(SM) # construction of independent reactions
3-element Vector{Reaction{Species{Int64}, Int64, Species{Int64}, Int64, Int64}}:
 H₂O = OH⁻ + H⁺
 H₂O + CO₂ = HCO₃⁻ + H⁺
 H₂O + CO₂ = CO₃²⁻ + 2H⁺
```

- Self-balancing of a chemical reaction: symbolic example of alkane combustion

```julia
julia> using SymPy

julia> n = symbols("n", real=true) ;

julia> CₙH₂ₙ₊₂ = Species(:C => n, :H => 2n+2) ;

julia> O₂, H₂O, CO₂ = Species.(split("O₂ H₂O CO₂")) ;

julia> r = Reaction([CₙH₂ₙ₊₂, O₂], [H₂O, CO₂])
CₙH₂ₙ₊₂ + (3n/2+1/2)O₂ = (n+1)H₂O + nCO₂
 reactants: CₙH₂ₙ₊₂ => 1, O₂ => 3*n/2 + 1/2
  products: H₂O => n + 1, CO₂ => n
    charge: 0

julia> pprint(r)
CₙH₂ₙ₊₂ + (3n/2+1/2)O₂ = (n+1)H₂O + nCO₂
 reactants: CₙH₂ₙ₊₂ => 1, O₂ => 3*n/2 + 1/2
  products: H₂O => n + 1, CO₂ => n
    charge: 0

julia> pprint(2r)
2CₙH₂ₙ₊₂ + (3n+1)O₂ = (2n+2)H₂O + 2nCO₂
 reactants: CₙH₂ₙ₊₂ => 2, O₂ => 3*n + 1
  products: H₂O => 2*n + 2, CO₂ => 2*n
    charge: 0

julia> for vn in 1:9 print("n=$vn ⇒ "); println(colored(apply(subs, r, n=>vn))) end
n=1 ⇒ CH₄ + 2O₂ = 2H₂O + CO₂
n=2 ⇒ C₂H₆ + 7/2O₂ = 3H₂O + 2CO₂
n=3 ⇒ C₃H₈ + 5O₂ = 4H₂O + 3CO₂
n=4 ⇒ C₄H₁₀ + 13/2O₂ = 5H₂O + 4CO₂
n=5 ⇒ C₅H₁₂ + 8O₂ = 6H₂O + 5CO₂
n=6 ⇒ C₆H₁₄ + 19/2O₂ = 7H₂O + 6CO₂
n=7 ⇒ C₇H₁₆ + 11O₂ = 8H₂O + 7CO₂
n=8 ⇒ C₈H₁₈ + 25/2O₂ = 9H₂O + 8CO₂
n=9 ⇒ C₉H₂₀ + 14O₂ = 10H₂O + 9CO₂
```

- Extraction from database ([ThermoHub]("https://github.com/thermohub") `.json` file) and reconstruction of thermodynamical functions

```julia
julia> json_file = "psinagra-12-07-thermofun.json"

julia> df_substances = read_thermofun_substances(json_file);

julia> dict_species = Dict(zip(df_substances.symbol, df_substances.species)); # dictionary indexing species by their symbol

julia> CaSO₄ = dict_species["Ca(SO4)@"]
Species{Int64}
           name: CaSO4  aq
         symbol: Ca(SO4)@
        formula: CaSO4@ ◆ CaSO₄@
          atoms: Ca => 1, S => 1, O => 4
         charge: 0
aggregate_state: AS_AQUEOUS
          class: SC_AQSOLUTE
     properties: M = 0.136133999952955 kg mol⁻¹
                 Tref = 298.15 K
                 Pref = 100000.0 m⁻¹ kg s⁻²
                 Cp⁰ = -104.60063934326 ♢ unit=[m² kg s⁻² K⁻¹ mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]
                 ΔfH⁰ = -1.417243319379807e6 - 104.60063934326T ♢ unit=[m² kg s⁻² mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]
                 S⁰ = 616.8922592448816 - 104.60063934326log(T) ♢ unit=[m² kg s⁻² K⁻¹ mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]
                 ΔfG⁰ = -1.27295402135706e6 - 721.4928985881417T + 104.60063934326T*log(T) ♢ unit=[m² kg s⁻² mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]
                 Vm = 4.700633883476301e-6 ♢ unit=[m³ mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]

julia> Ca²⁺ = dict_species["Ca+2"]
Species{Int64}
           name: Ca+2
         symbol: Ca+2
        formula: Ca+2 ◆ Ca²⁺
          atoms: Ca => 1
         charge: 2
aggregate_state: AS_AQUEOUS
          class: SC_AQSOLUTE
     properties: M = 0.04007799998614991 kg mol⁻¹
                 Tref = 298.15 K
                 Pref = 100000.0 m⁻¹ kg s⁻²
                 Cp⁰ = -30.922515869141 ♢ unit=[m² kg s⁻² K⁻¹ mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]
                 ΔfH⁰ = -533849.4518936157 - 30.922515869141T ♢ unit=[m² kg s⁻² mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]
                 S⁰ = 119.70002369348362 - 30.922515869141log(T) ♢ unit=[m² kg s⁻² K⁻¹ mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]
                 ΔfG⁰ = -560411.1568393706 - 150.6225395626246T + 30.922515869141T*log(T) ♢ unit=[m² kg s⁻² mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]
                 Vm = -1.8438742160797e-5 ♢ unit=[m³ mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]

julia> SO₄²⁻ = dict_species["SO4-2"]
Species{Int64}
           name: SO4-2
         symbol: SO4-2
        formula: S|6|O4-2 ◆ SO₄²⁻
          atoms: S => 1, O => 4
         charge: -2
aggregate_state: AS_AQUEOUS
          class: SC_AQSOLUTE
     properties: M = 0.09605599996680511 kg mol⁻¹
                 Tref = 298.15 K
                 Pref = 100000.0 m⁻¹ kg s⁻²
                 Cp⁰ = -266.09072875977 ♢ unit=[m² kg s⁻² K⁻¹ mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]
                 ΔfH⁰ = -830362.0492202746 - 266.09072875977T ♢ unit=[m² kg s⁻² mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]
                 S⁰ = 1534.9056613400478 - 266.09072875977log(T) ♢ unit=[m² kg s⁻² K⁻¹ mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]
                 ΔfG⁰ = -659510.4812841403 - 1800.9963900998177T + 266.09072875977T*log(T) ♢ unit=[m² kg s⁻² mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]
                 Vm = 1.2917655706406e-5 ♢ unit=[m³ mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]

julia> r = Reaction([CaSO₄, Ca²⁺, SO₄²⁻])
  equation: CaSO₄@ = Ca²⁺ + SO₄²⁻
 reactants: CaSO₄@ => 1
  products: Ca²⁺ => 1, SO₄²⁻ => 1
    charge: 0
properties: ΔrCp⁰ = -192.412605285651 ♢ unit=[m² kg s⁻² K⁻¹ mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]
            ΔrS⁰ = 1037.7134257886498 - 192.412605285651log(T) ♢ unit=[m² kg s⁻² K⁻¹ mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]
            ΔrH⁰ = 53031.81826591678 - 192.412605285651T ♢ unit=[m² kg s⁻² mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]
            ΔrG⁰ = 53032.38323354907 - 1230.1260310743007T + 192.412605285651T*log(T) ♢ unit=[m² kg s⁻² mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]
            ΔrV = -1.02217203378673e-5 ♢ unit=[m³ mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]
```

## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add ChemistryLab
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("ChemistryLab")
```

## Usage

See the [documentation and tutorials](https://jfbarthelemy.github.io/ChemistryLab.jl) for examples on formula creation, species management, reaction parsing, and database merging.

## License

MIT License. See [LICENSE](LICENSE) for details.

## Citation

[![DOI](https://zenodo.org/badge/1054296488.svg)](https://doi.org/10.5281/zenodo.17756074)

See [CITATION.cff](CITATION.cff) for citation details.

**BibTeX entry:**

```bibtex
@software{chemistrylab_jl,
  authors = {Barthélemy, Jean-François and Soive, Anthony},
  title = {ChemistryLab.jl: Numerical laboratory for computational chemistry},
  doi = {10.5281/zenodo.17756074},
  url = {https://github.com/jfbarthelemy/ChemistryLab.jl}
}
```

## Credits and Acknowledgements

Developed by [Jean-François Barthélémy](https://github.com/jfbarthelemy) and [Anthony Soive](https://github.com/anthonysoive), both researchers at [Cerema](https://www.cerema.fr/en) in the research team [UMR MCD](https://mcd.univ-gustave-eiffel.fr/).
