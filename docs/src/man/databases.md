# Database Interoperability

So far, we have looked at the possibility of creating and manipulating any species, whether they exist or not. If we wanted to create a H₂O⁺⁴ molecule, it would not be a problem. However, you will admit that it is a little strange...

## Cemdata18 and PSI-Nagra-12-07 Databases

This is why Cement Chemistry relies on existing databases, in particular [Cemdata18](https://www.empa.ch/web/s308/thermodynamic-data) and [PSI-Nagra-12-07](https://www.psi.ch/en/les/thermodynamic-databases). Cemdata18 is a chemical thermodynamic database for hydrated Portland cements and alkali-activated materials. PSI-Nagra is a Chemical Thermodynamic Database. The formalism adopted for these databases is that of [Thermofun](https://thermohub.org/thermofun/thermofun/) which is a universal open-source client that delivers thermodynamic properties of substances and reactions at the temperature and pressure of interest. The information is stored in json files.

With ChemistryLab, you can parse a ThermoFun-like json file and return DataFrames for:

- elements:

```@example database_interoperability
using ChemistryLab #hide
df_elements, df_substances, df_reactions = read_thermofun("../../../data/cemdata18-merged.json")
show(df_elements, allcols=true, allrows=true)
```
[`ChemistryLab.read_thermofun`](@ref)

- species (aqueous, solid or gaseous phases):

```@example database_interoperability
show(df_substances, allcols=true, allrows=false)
```

- reactions contained in the database
```@example database_interoperability
show(df_reactions, allcols=true, allrows=false)
```

!!! note "Units"
    By default, when reading the database, the units of the various physical quantities are read. To remove this option, simply write `with_units=false` in the arguments of the `read_thermofun` function.

## Species properties construction during database parsing
Species properties can be constructed while reading the database with `all_properties=true` as an argument. These properties are values ​​or functions that generally depend on temperature.

Currently, species properties are: molar mass (`molar_mass`), a reference temperature (`Tref`), heat capacity (`Cp`), enthalpy change of formation (`ΔfH`), entropy of formation (`S`), Gibbs energy of formation (`ΔfG`), and molar volume (`Vm`). Values and functions are described in [Cemdata18 paper](https://www.empa.ch/web/s308/thermodynamic-data).

```@example database_interoperability
df_elements, df_substances, df_reactions = read_thermofun("../../../data/cemdata18-merged.json"; with_units=true, all_properties=true) #hide
```

For each species, properties can then be obtained as follows:
```@example database_interoperability
dict_species = Dict(zip(df_substances.symbol, df_substances.species))
dict_species["Portlandite"]
```

## Primary species extraction

It is also possible to retrieve primary species from the Cemdata18 database, primary species being the designation of a subset of species for which any species can be represented as the linear combination of primary species.

```@example database_interoperability
df_primaries = extract_primary_species("../../../data/CEMDATA18-31-03-2022-phaseVol.dat")
show(df_primaries, allcols=true, allrows=true)
```
[`ChemistryLab.extract_primary_species`](@ref)

---

