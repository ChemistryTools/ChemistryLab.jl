<!-- ![ChemistryLab Logo](https://raw.githubusercontent.com/jfbarthelemy/ChemistryLab.jl/main/docs/src/assets/logo.svg) -->
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

## Installation

```julia
using Pkg
Pkg.add("ChemistryLab")
```

## Usage

See the [documentation and tutorials](https://jfbarthelemy.github.io/ChemistryLab.jl) for examples on formula creation, species management, reaction parsing, and database merging.

## License

MIT License. See [LICENSE](LICENSE) for details.

## Credits and Acknowledgements

Developed by [Jean-François Barthélémy](https://github.com/jfbarthelemy) and [Anthony Soive](https://github.com/anthonysoive), both researchers at [Cerema](https://www.cerema.fr/en) in the research team [UMR MCD](https://mcd.univ-gustave-eiffel.fr/).
