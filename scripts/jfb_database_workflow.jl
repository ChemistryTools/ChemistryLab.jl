using ChemistryLab

df_elements, df_substances, df_reactions = read_thermofun_database("data/cemdata18-merged.json")
df_union = get_compatible_species(split("C3S C2S C3A Portlandite Jennite H2O@"), df_substances;
               aggregate_states=[AS_AQUEOUS], exclude_species=split("H2@ O2@"), union=true)
dict_all_species = build_species_from_database(df_union)
candidate_primaries = [dict_all_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_all_species, s)]
SM = StoichMatrix(dict_all_species, candidate_primaries) ; pprint(SM)
list_reactions = reactions(SM)
pprint(list_reactions)
for r in list_reactions
    display(r)
end
