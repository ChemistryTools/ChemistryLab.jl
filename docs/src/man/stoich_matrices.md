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

In chemistryLab, it is possible to construct a stoichiometric matrix of species as a function of primary species given a database. Primary species candidates can be found in a database. Those from Cemdata18 can be listed with the following command:

```julia
using ChemistryLab
df_elements, df_substances, df_reactions = read_thermofun("../../../data/cemdata18-merged.json")
df_primaries = extract_primary_species("../../../data/CEMDATA18-31-03-2022-phaseVol.dat")
dict_species = Dict(zip(df_substances.symbol, df_substances.species))
candidate_primaries = [s == "Zz" ? Species("Zz") : dict_species[s] for s in df_primaries.symbol]
```

For a given number of species (eg. Portlandite and water), all ionic species which could appear during chemical reactions have to be identified:

```julia
given_species = filter(row -> row.symbol ∈ split("Portlandite H2O@"), df_substances)
secondaries = filter(row->row.aggregate_state == "AS_AQUEOUS" 
                          && all(k->first(k) ∈ union_atoms(atoms.(given_species.species)), atoms(row.species))
                          && row.symbol ∉ split("H2@ O2@"),
                          df_substances)
```

Primary species concerning by the reactions can then be deduced.

```julia
species = unique(vcat(given_species, secondaries), :symbol).species
candidate_primaries = [s == "Zz" ? Species("Zz") : dict_species[s] for s in df_primaries.symbol]
```

Finally, the stoichiometric matrix can be calculated:

```@setup example1
    using ChemistryLab #hide
    # using Serialization
    using PrettyTables
    # df_substances, df_reactions = deserialize("../../../data/cemdata18.jls")
    df_elements, df_substances, df_reactions = read_thermofun("../../../data/cemdata18-merged.json") #hide
    df_primaries = extract_primary_species("../../../data/CEMDATA18-31-03-2022-phaseVol.dat") #hide
    dict_species = Dict(zip(df_substances.symbol, df_substances.species)) #hide

    given_species = filter(row -> row.symbol ∈ split("Portlandite H2O@"), df_substances) #hide
    secondaries = filter(row->row.aggregate_state == "AS_AQUEOUS" 
                            && all(k->first(k) ∈ union_atoms(atoms.(given_species.species)), atoms(row.species))
                            && row.symbol ∉ split("H2@ O2@"),
                            df_substances)


    species = unique(vcat(given_species, secondaries), :symbol).species #hide
    candidate_primaries = [s == "Zz" ? Species("Zz") : dict_species[s] for s in df_primaries.symbol] #hide
```

```@example example1
SM = StoichMatrix(species, candidate_primaries)

using PrettyTables #hide
```

---
