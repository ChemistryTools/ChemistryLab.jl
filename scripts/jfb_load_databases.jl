using Revise
using ChemistryLab, Unicode
using DynamicQuantities
using ModelingToolkit
using JSON
using DataFrames

# println("LOADING DATABASES...")
# @time df_elements, df_substances, df_reactions, dict_species, dict_reactions = read_thermofun("data/cemdata18-merged")
# @time df_elements_psi, df_substances_psi, df_reactions_psi, dict_species_psi, dict_reactions_psi = read_thermofun("data/psinagra-12-07-thermofun")
# @time df_elements_aq17, df_substances_aq17, df_reactions_aq17, dict_species_aq17, dict_reactions_aq17 = read_thermofun("data/aq17-thermofun")
# @time df_elements_orga, df_substances_orga, df_reactions_orga, dict_species_orga, dict_reactions_orga = read_thermofun("data/slop98-organic-thermofun")

filename = "data/cemdata18-merged.json"
json_str = read(filename, String)
data = JSON.parse(json_str)
data = JSON.parsefile(filename)
df_substances = DataFrame(Tables.dictrowtable(data["substances"]))
df_substances.TPMethods[10][1]["eos_hkf_coeffs"]["units"]

df_substances = DataFrame(Tables.dictrowtable(data["substances"]))
df_reactions = DataFrame(Tables.dictrowtable(data["reactions"]))
df_elements = DataFrame(Tables.dictrowtable(data["elements"]))
list_symbols = split("C2S Portlandite Jennite H2O@")
local_df_substances = filter(row -> row.symbol ∈ list_symbols, df_substances)
row = eachrow(local_df_substances)[4]
TPMethods = row.TPMethods
method = TPMethods[1]
method_type = only(values(method.method))
coeffs = method.m_heat_capacity_ft_coeffs
uparse.(coeffs.units)
coeffs.values

row = eachrow(df_reactions)[4]

println("LOADING DATABASES...")
df_elements, df_substances, df_reactions = read_thermofun_database("data/cemdata18-merged.json")
species_vec = build_species(df_substances)
dict_species = Dict(symbol(s) => s for s in species_vec)
dict_reactions = Dict(symbol(r) => r for r in build_reactions(df_reactions, species_vec))

df_elements_psi, df_substances_psi, df_reactions_psi = read_thermofun_database("data/psinagra-12-07-thermofun.json")
species_psi = build_species(df_substances_psi)
dict_species_psi = Dict(symbol(s) => s for s in species_psi)
dict_reactions_psi = Dict(symbol(r) => r for r in build_reactions(df_reactions_psi, species_psi))

df_elements_aq17, df_substances_aq17, df_reactions_aq17 = read_thermofun_database("data/aq17-thermofun.json")
species_aq17 = build_species(df_substances_aq17)
dict_species_aq17 = Dict(symbol(s) => s for s in species_aq17)
dict_reactions_aq17 = Dict(symbol(r) => r for r in build_reactions(df_reactions_aq17, species_aq17))

df_elements_orga, df_substances_orga, df_reactions_orga = read_thermofun_database("data/slop98-organic-thermofun.json")
species_orga = build_species(df_substances_orga)
dict_species_orga = Dict(symbol(s) => s for s in species_orga)
dict_reactions_orga = Dict(symbol(r) => r for r in build_reactions(df_reactions_orga, species_orga))

df_elements_inorga, df_substances_inorga, df_reactions_inorga = read_thermofun_database("data/slop98-inorganic-thermofun.json")
species_inorga = build_species(df_substances_inorga)
dict_species_inorga = Dict(symbol(s) => s for s in species_inorga)
dict_reactions_inorga = Dict(symbol(r) => r for r in build_reactions(df_reactions_inorga, species_inorga))
