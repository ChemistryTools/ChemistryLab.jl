# Database Interoperability

So far, we have looked at the possibility of creating and manipulating any species, whether they exist or not. If we wanted to create a H₂O⁺⁴ molecule, it would not be a problem. However, you will admit that it is a little strange...

This is why ChemistryLab relies on existing databases, in particular [Cemdata18](https://www.empa.ch/web/s308/thermodynamic-data) and [PSI-Nagra-12-07](https://www.psi.ch/en/les/thermodynamic-databases). Cemdata18 is a chemical thermodynamic database for hydrated Portland cements and alkali-activated materials. PSI-Nagra is a Chemical Thermodynamic Database. The formalism adopted for these databases is that of [Thermofun](https://thermohub.org/thermofun/thermofun/) which is a universal open-source client that delivers thermodynamic properties of substances and reactions at the temperature and pressure of interest. The information is stored in json files.

## Loading species from a database

The simplest way to load species from a ThermoFun-compatible JSON file is `build_species`, which reads the file and directly returns a `Vector{Species}` with compiled thermodynamic functions:

```julia
using ChemistryLab
all_species = build_species("../../../data/cemdata18-merged.json")
```

Each species already carries its molar mass and temperature-dependent thermodynamic functions (Cp⁰, ΔₐH⁰, ΔₐS⁰, ΔₐG⁰, logK⁰) as `ThermoFunction`s.

!!! note "Low-level access"
    If you need the raw DataFrames (e.g. to inspect metadata or filter on database columns), the lower-level function `read_thermofun_database` is still available and returns three DataFrames `(df_elements, df_substances, df_reactions)`:
    ```julia
    df_elements, df_substances, df_reactions = read_thermofun_database("../../../data/cemdata18-merged.json")
    ```
    `build_species(df_substances)` can then be called on the filtered DataFrame.

## Filtering species with `speciation`

In practice, only a small subset of the database is relevant to a given problem. `speciation` filters a species list to those whose atomic composition is a subset of the atoms found in a set of *seed* species:

```@example database
using ChemistryLab #hide
all_species = build_species("../../../data/cemdata18-merged.json") #hide
# Keep only species that can form from the calcite / water system
species_calcite = speciation(all_species, split("Cal H2O@");
                             aggregate_state=[AS_AQUEOUS],
                             exclude_species=split("H2@ O2@ CH4@"))
dict_species_calcite = Dict(symbol(s) => s for s in species_calcite)
```

The `aggregate_state` keyword restricts results to aqueous species (the seed species themselves — `Cal` and `H2O@` — are always kept through `include_species` internally). For example, the properties of Ca(HCO₃)⁺ can then be read as:

```@example database
dict_species_calcite["Ca(HCO3)+"]
```

### `speciation` signatures

`speciation` accepts seed arguments in three forms:

| Seed argument | Description |
|---------------|-------------|
| `Vector{Symbol}` | Explicit list of atom symbols |
| `Vector{<:AbstractSpecies}` | Species objects — their union of atoms defines the space |
| `Vector{<:AbstractString}` | Species symbol strings — looked up in `species_list` |

Common keyword arguments:

| Keyword | Default | Description |
|---------|---------|-------------|
| `aggregate_state` | all states | restrict to `[AS_AQUEOUS]`, `[AS_CRYSTAL]`, etc. |
| `class` | all classes | restrict to `[SC_AQSOLUTE]`, etc. |
| `exclude_species` | `[]` | species (or symbols) to always exclude |
| `include_species` | `[]` | species to always include regardless of composition |

## Primary species extraction

It is also possible to retrieve primary species from the Cemdata18 database. Primary species are a minimal subset such that every other species can be expressed as their linear combination.

```julia
df_primaries = extract_primary_species("../../../data/CEMDATA18-31-03-2022-phaseVol.dat")
show(df_primaries, allcols=true, allrows=true)
```

---
