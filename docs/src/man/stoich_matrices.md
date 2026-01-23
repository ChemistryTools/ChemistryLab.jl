# Stoichiometric Matrix

From the definition of species, it is possible to construct a stoichiometric matrix that establishes the relationship between species and chemical elements for species or oxides for cement species. This is called canonical decomposition.

```@setup database_stoichiometry
    using ChemistryLab
    import Pkg; Pkg.add("PrettyTables")
```

## Stoichiometric matrix for species

Any species can be described as a linear combination of chemical elements. A species vector can be expressed as a function of the chemical elements on which they depend. This dependence leads to the creation of a stoichiometric matrix.

```@example
using ChemistryLab #hide
H2O = Species("H₂O")
HSO4 = Species("HSO₄⁻")
CO2 = Species(Dict(:C => 1, :O => 2); symbol="CO₂")
species = [H2O, HSO4, CO2]
CSM = CanonicalStoichMatrix(species)

using PrettyTables #hide
```

## Stoichiometric matrix for cement species

A cement species vector can also be expressed in terms of other species on which they depend. Here, the cement species are expressed in terms of the oxides from which they are composed.

```@example stoich
using ChemistryLab #hide
C3S = CemSpecies("C3S")
C2S = CemSpecies("C2S")
C3A = CemSpecies("C3A")
C4AF = CemSpecies(Dict(:C=>4, :A=>1, :F=>1); name="C4AF")
cemspecies = [C3S, C2S, C3A, C4AF]
CSM = CanonicalStoichMatrix(cemspecies)

using PrettyTables #hide
```

---

## Define stoichiometric matrix from primary species

The decomposition of a set of species can also be done according to a base of primary species.

```@example stoich
using ChemistryLab #hide
H2O = Species("H₂O")
HSO4 = Species("HSO₄⁻")
CO2 = Species(Dict(:C => 1, :O => 2); symbol="CO₂")
species = [H2O, HSO4, CO2]
H⁺ = Species("H⁺")
SO₄²⁻ = Species("SO₄²⁻")
CO₃²⁻ = Species("CO₃²⁻")
primary_species = [H⁺, SO₄²⁻, CO₃²⁻, H2O]
SM = StoichMatrix(species, primary_species)

using PrettyTables #hide
```

!!! note "Display of the stoichiometric matrix"
    The stoichiometric matrix can be pretty-printed with different column and row labels using `pprint`. Simply add the keyword `row_label`, `col_label` or `label` for both, which can take the following values: *:name*, *:symbol*, *:formula*
    ```julia
    SM = StoichMatrix(species)
    pprint(SM; label=:name)
    ```

## Construct stoichiometric matrix from database

In chemistryLab, it is possible to construct a stoichiometric matrix of species as a function of primary species given a database. First the database is loaded by

```julia
using ChemistryLab
df_elements, df_substances, df_reactions = read_thermofun_database("../../../data/cemdata18-merged.json")
```

For a given number of species (eg. Portlandite and water), all ionic species which could appear during chemical reactions have to be identified. Here all species of the database entirely composed of atoms involved in the given species under provided aggregate states and possibly excluding a selection of species are gathered in a `DataFrame` obtained from `get_compatible_species` as

```julia
df_given_species = filter(row -> row.symbol ∈ split("Portlandite H2O@"), df_substances)
df_added_species = get_compatible_species(split("Portlandite H2O@"), df_substances;
               aggregate_states=[AS_AQUEOUS], exclude_species=split("H2@ O2@"))
```

The given and added species can be grouped in a unique `DataFrame` by

```julia
df_union = unique(vcat(df_given_species, df_added_species))
```

Or in one line with the keyword argument `union=true`

```julia
df_union = get_compatible_species(split("Portlandite H2O@"), df_substances;
               aggregate_states=[AS_AQUEOUS], exclude_species=split("H2@ O2@ FeOH+ Fe+2"), union=true)
```

A dictionary of `Species` objects can then be built from a `DataFrame` thanks to:

```julia
dict_species = build_species_from_database(df_union)
```

!!! warning "Calculation of Molar Mass"
    The function `build_species_from_database` is rather time costly due to the construction of `ThermoFunction`s so that it can be long on large databases.

Primary species candidates can be found in a database. Those from Cemdata18 can be listed with the following command:

```julia
candidate_primaries = [dict_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_species, s)]
```

Finally, the stoichiometric matrix can be calculated:

```@setup example1
    using ChemistryLab #hide
    df_elements, df_substances, df_reactions = read_thermofun_database("../../../data/cemdata18-merged.json") #hide
    df_union = get_compatible_species(split("Portlandite H2O@"), df_substances;
                aggregate_states=[AS_AQUEOUS], exclude_species=split("H2@ O2@ FeOH+ Fe+2"), union=true) #hide
    dict_species = build_species_from_database(df_union) #hide
    candidate_primaries = [dict_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_species, s)] #hide
```

```@example example1
SM = StoichMatrix(dict_species, candidate_primaries)
```

---
