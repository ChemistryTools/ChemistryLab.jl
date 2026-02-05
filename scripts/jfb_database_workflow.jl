using Revise
using ChemistryLab
using OrderedCollections

df_elements, df_substances, df_reactions = read_thermofun_database("data/cemdata18-merged.json")
df_union = get_compatible_species(split("C3S C2S C3A C4AF Portlandite Jennite H2O@"), df_substances;
               aggregate_states=[AS_AQUEOUS], exclude_species=split("H2@ O2@ FeOH+ Fe+2"), union=true)
dict_all_species = build_species_from_database(df_union)
candidate_primaries = [dict_all_species[s] for s in CEMDATA_PRIMARIES if haskey(dict_all_species, s)]
SM = StoichMatrix(dict_all_species, candidate_primaries) ; pprint(SM)
list_reactions = reactions(SM) ; pprint(list_reactions)
for r in list_reactions display(r); println() end

df_hydrates_clinker = get_compatible_species(split("C3S C2S C3A C4AF Gp H2O@"), df_substances; aggregate_states=[AS_CRYSTAL, AS_AQUEOUS])
dict_hydrates = build_species_from_database(df_hydrates_clinker)
candidate_primaries = [dict_hydrates[s] for s in CEMDATA_PRIMARIES if haskey(dict_hydrates, s)]
@time SM = StoichMatrix(dict_hydrates, candidate_primaries) ; pprint(SM)
@time list_reactions = reactions(SM) ; pprint(list_reactions)

dict_species = build_species_from_database(df_substances)
dict_reactions = build_reactions_from_database(df_reactions, dict_species)

F = eval(build_function(func, [T]; expression=Val{false}))
Cp = Num(parse_expr_to_symbolic(THERMO_MODELS[:cp_ft_equation][:Cp], @__MODULE__))

params = [:a₀ => 210.0, :a₁ => 0.12, :a₂ => -3.07e6, :a₃ => 0.0, :T => 298.15, :Cp⁰ => 210.0, :ΔfH⁰ => -2723484.33, :S⁰ => 140, :ΔfG⁰ => -2480808.197, :V⁰ => 7.84]
dict_params = Dict(Num(parse_expr_to_symbolic(k, @__MODULE__)) => v for (k,v) in params)
dict_params = Dict(k => get(dict_params, k, 0) for k in get_variables(Cp))
ModelingToolkit.substitute(Cp, dict_params)
