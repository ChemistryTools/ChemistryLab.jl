using Revise, ChemistryLab, Unicode
using DynamicQuantities
using ModelingToolkit
using Serialization
import SymPy: symbols, Sym, N, subs, factor, simplify

# json_file, jls_file = "data/cemdata18-merged.json", "data/cemdata18.jls"
json_file, jls_file = "data/psinagra-12-07-thermofun.json", "data/psinagra.jls"
try
    global df_substances, df_reactions = deserialize(jls_file)
catch
    global df_substances = read_thermofun_substances(json_file; with_units=true, add_species=true, all_properties=true, debug=false)
    global df_reactions = read_thermofun_reactions(json_file, df_substances; with_units=true, add_reactions=true, all_properties=true, debug=false)
    serialize(jls_file, (df_substances, df_reactions))
end
dict_species = Dict(zip(df_substances.symbol, df_substances.species))
dict_reactions = Dict(zip(df_reactions.symbol, df_reactions.reaction))
