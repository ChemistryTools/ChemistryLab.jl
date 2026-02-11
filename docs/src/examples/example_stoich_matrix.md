# Stoichiometric Matrix

Calculating stoichiometric matrices is a prerequisite for equilibrium calculations by minimizing Gibbs energy. The examples below show how they can be constructed.

## Get Stoichiometric Matrix from a list of species

Let's imagine that we want to form the stoichiometric matrix of a list of solid and water species. For that, we need to read the database from which these species originate.

```julia
using ChemistryLab
using PrettyTables
df_elements, df_substances, df_reactions = read_thermofun_database("../../../data/cemdata18-merged.json")
```

It is then necessary to identify the list of secondary species likely to appear during the reactions.

```julia
df_species = get_compatible_species(df_substances, split("C3S Portlandite Jennite H2O@");
               aggregate_states=[AS_AQUEOUS], exclude_species=split("H2@ O2@"), union=true)
dict_species = build_species_from_database(df_species)
```

We can then deduce the primary species concerned by the reaction.

```julia
candidate_primaries = [dict_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_species, s)]
```

And construct the stoichiometric matrix

```@setup example1
    using ChemistryLab #hide
    df_elements, df_substances, df_reactions = read_thermofun_database("../../../data/cemdata18-merged.json") #hide

    df_species = get_compatible_species(df_substances, split("C3S Portlandite Jennite H2O@");
                aggregate_states=[AS_AQUEOUS], exclude_species=split("H2@ O2@"), union=true) #hide
    dict_species = build_species_from_database(df_species) #hide

    candidate_primaries = [dict_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_species, s)] #hide
```

```@example example1
SM = StoichMatrix(dict_species, candidate_primaries)
```

## Get Stoichiometric Matrix from a database file

```julia
using ChemistryLab

df_elements, df_substances, df_reactions = read_thermofun_database("../../../data/cemdata18-merged.json")
```

```@example example1
df_aqueous_species = filter(row -> only(row.aggregate_state).second == "AS_AQUEOUS", df_substances)
dict_aqueous_species = build_species_from_database(df_aqueous_species)
candidate_primaries = [dict_aqueous_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_aqueous_species, s)]
SM = StoichMatrix(dict_aqueous_species, candidate_primaries)
```

All the independent reactions of the species contained in the database can thus be reconstructed. Here, only ionic species are listed given the choice to only read ionic species in the database ("AS_AQUEOUS").

```@example example1
lr = reactions(SM)
```

---

The exercise can also be done on solid species. In this case, the data filter is carried out using the keyword "AS_CRYSTAL", in accordance with the terminology adopted in Thermofun. Note that primaries are still built among aqueous species.

```@setup example1
df_solid_species = filter(row -> only(row.aggregate_state).second == "AS_CRYSTAL", df_substances)
dict_solid_species = build_species_from_database(df_solid_species)
candidate_primaries = [dict_aqueous_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_aqueous_species, s)]
SM = StoichMatrix(dict_solid_species, candidate_primaries) ; pprint(SM)
```

```@example example1
lr = reactions(SM)
```

---

Or with gases ("AS_GAS")

```@setup example1
df_gas_species = filter(row -> only(row.aggregate_state).second == "AS_GAS", df_substances)
dict_gas_species = build_species_from_database(df_gas_species)
candidate_primaries = [dict_aqueous_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_aqueous_species, s)]
SM = StoichMatrix(dict_gas_species, candidate_primaries) ; pprint(SM)
```

```@example example1
lr = reactions(SM)
```
