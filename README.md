<p>
  <img src="./docs/src/assets/logo.svg" width="100">
</p>

# ChemistryLab

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jfbarthelemy.github.io/ChemistryLab.jl/dev/)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jfbarthelemy.github.io/ChemistryLab.jl/stable/)
[![Build Status](https://github.com/jfbarthelemy/ChemistryLab.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jfbarthelemy/ChemistryLab.jl/actions/workflows/CI.yml?query=branch%3Amain)

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

## Credits and Acknowledgements

Developed by [Jean-François Barthélémy](https://github.com/jfbarthelemy) and [Anthony Soive](https://github.com/anthonysoive), both researchers at [Cerema](https://www.cerema.fr/en) in the research team [UMR MCD](https://mcd.univ-gustave-eiffel.fr/).
