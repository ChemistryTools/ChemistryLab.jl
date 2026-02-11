# Database Interoperability

So far, we have looked at the possibility of creating and manipulating any species, whether they exist or not. If we wanted to create a H₂O⁺⁴ molecule, it would not be a problem. However, you will admit that it is a little strange...

This is why ChemistryLab relies on existing databases, in particular [Cemdata18](https://www.empa.ch/web/s308/thermodynamic-data) and [PSI-Nagra-12-07](https://www.psi.ch/en/les/thermodynamic-databases). Cemdata18 is a chemical thermodynamic database for hydrated Portland cements and alkali-activated materials. PSI-Nagra is a Chemical Thermodynamic Database. The formalism adopted for these databases is that of [Thermofun](https://thermohub.org/thermofun/thermofun/) which is a universal open-source client that delivers thermodynamic properties of substances and reactions at the temperature and pressure of interest. The information is stored in json files.

## Database parsing

With ChemistryLab, you can parse a ThermoFun-like json file and return DataFrames for elements, species (aqueous, solid or gaseous phases) or reactions:

```julia
using ChemistryLab
filebasename = "cemdata18-merged.json"
df_elements, df_substances, df_reactions = read_thermofun_database("../../../data/" * filebasename)
```

For easier reading, the following commands can be run:

```@example database
using ChemistryLab #hide
filebasename = "cemdata18-merged.json" #hide
df_elements, df_substances, df_reactions = read_thermofun_database("../../../data/" * filebasename) #hide
show(df_elements, allcols=true, allrows=true)
show(df_substances, allcols=true, allrows=false)
show(df_reactions, allcols=true, allrows=false)
```

## Species properties construction during database parsing
All species in the database could be constructed from Dataframes obtained by reading the database. However, this operation is time-consuming and of little practical value, as the chemical reactions of interest are often a much smaller subset.

Let's take a subset representing the dissolution of calcite in pure water. `get_compatible_species` allows us to construct the subset of species that can be present during the reaction of calcite in water.

```julia
df_calcite = get_compatible_species(split("Cal H2O@"), df_substances;
                        aggregate_states=[AS_AQUEOUS], exclude_species=split("H2@ O2@ CH4@"), union=true)
```

The construction of thermodynamic functions is then done by calling the `build_species_from_database` function:

```julia
dict_species_calcite = build_species_from_database(df_calcite)
```

For example, $\text{Ca}(\text{HCO}_3)^+$ properties can be read as follows:

```@example database
df_calcite = get_compatible_species(split("Cal H2O@"), df_substances;
                        aggregate_states=[AS_AQUEOUS], exclude_species=split("H2@ O2@ CH4@"), union=true) #hide
dict_species_calcite = build_species_from_database(df_calcite) #hide
dict_species_calcite["Ca(HCO3)+"]
```

## Primary species extraction

It is also possible to retrieve primary species from the Cemdata18 database, primary species being the designation of a subset of species for which any species can be represented as the linear combination of primary species.

```julia
df_primaries = extract_primary_species("../../../data/CEMDATA18-31-03-2022-phaseVol.dat")
show(df_primaries, allcols=true, allrows=true)
```

---
