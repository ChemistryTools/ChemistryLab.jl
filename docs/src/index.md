```@meta
CurrentModule = ChemistryLab
```

# ChemistryLab

Welcome to the ChemistryLab documentation. This set of pages guides you through the core features of ChemistryLab — a Julia toolkit for parsing and manipulating chemical formulas, creating species, building stoichiometric matrices and reactions, and importing thermodynamic data (ThermoFun / Cemdata).

## Audience

This documentation is intended for scientists, engineers and scientific programmers who want programmatic control over chemical species and simple thermodynamic workflows in Julia. It suits chemists, geochemists, civil-materials researchers (cement chemistry), and anyone who needs precise, scriptable handling of formulas, reactions and databases.

## Installation

To install ChemistryLab.jl, use the Julia package manager:

- From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add ChemistryLab
```

- Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("ChemistryLab")
```

## Citation

If you use ChemistryLab in your work, please cite the following:

```
@software{barthelemy_2026_18216605,
  author       = {Barthélemy, Jean-François and
                  Soive, Anthony},
  title        = {ChemistryLab.jl: Numerical laboratory for
                   computational chemistry
                  },
  month        = jan,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v0.1.5},
  doi          = {10.5281/zenodo.18216605},
  url          = {https://doi.org/10.5281/zenodo.18216605},
}
```

## Feature

- Parse and manipulate chemical formulas (supporting charges, subscripts/superscripts and rational stoichiometry).
- Create and inspect `Species` and `CemSpecies` objects and attach thermodynamic or auxiliary properties.
- Build canonical stoichiometric matrices and convert them to `Reaction` objects.
- Read, merge and write ThermoFun / Cemdata sources and extract species and reactions programmatically.
- Combine, simplify and transform reactions using the `Reaction` API.

## Documentation structure

- `introduction` — this page: goals, audience, and how to run the examples.
- `quickstart` — minimal examples to get you started quickly.
- `tutorial` — basic examples to get you started more deeply in the library.
  - `formula_manipulation` — parsing, conversions (Phreeqc ⇄ Unicode), formatting and arithmetic on `Formula` objects.
  - `species` — constructing `Species` and `CemSpecies`, managing properties and basic queries.
  - `equations` — parsing, formatting and coloring of chemical equations and balancing helpers.
  - `stoich_matrices` — building canonical stoichiometric matrices and converting them to reaction objects.
  - `databases` — reading ThermoFun JSON and Cemdata `.dat` files, merging reactions and exporting results.
- `example/` — runnable examples and small worked problems (also built to HTML under `docs/build/example/`).

## Try it locally

The repository contains a `docs` folder with Documenter sources and a `make.jl` script. To build the documentation locally from the project root, activate the docs environment and run the make script:

```powershell
# from the project root (Windows PowerShell)
julia --project=docs docs/make.jl
```

Alternatively, ensure the package environment is prepared and run Documenter through the package project:

```powershell
julia --project=. -e "using Pkg; Pkg.instantiate(); include(\"docs/make.jl\")"
```

After the build completes, open `docs/build/index.html` in a browser to read the rendered tutorial and follow the example pages.

## Quick tips

- In the REPL try small calls like `using ChemistryLab; Species("CaCO3")` and `Formula("SO4-2")` to explore parsing behavior interactively.
- Start with `docs/src/example/get_stoichio_matrix.md` to see a concise, runnable example converting a stoichiometric matrix into reactions.
- If you plan to work with ThermoFun/Cemdata sources, run the examples in `docs/src/man/databases.md` after placing the required `.json`/`.dat` data files in the `data/` directory.

## Next steps

You can see `examples` for more advanced runnable examples and small worked problems.

Happy exploring — this tutorial aims to be practical and runnable, so please tell me which example you want expanded into a fully reproducible script.
