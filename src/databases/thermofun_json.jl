"""
    read_thermofun_database(filename::AbstractString) -> (DataFrame, DataFrame, DataFrame)

Read a ThermoFun database from a JSON file.

# Arguments

  - `filename`: path to the JSON database file.

# Returns

  - `df_elements`: DataFrame of chemical elements.
  - `df_substances`: DataFrame of chemical substances (species).
  - `df_reactions`: DataFrame of chemical reactions.
"""
function read_thermofun_database(filename)
    print_title(
        "Loading database: $filename";
        crayon=Crayon(; foreground=:green),
        style=:box,
        indent="",
    )
    data = JSON.parsefile(filename)
    df_substances = DataFrame(Tables.dictrowtable(data["substances"]))
    df_reactions = DataFrame(Tables.dictrowtable(data["reactions"]))
    df_elements = DataFrame(Tables.dictrowtable(data["elements"]))
    return df_elements, df_substances, df_reactions
end

function extract_unit(v, default_unit=u"1")
    return try
        uparse(v)
    catch
        default_unit
    end
end

function extract_value(
    row, field::Symbol; verbose=false, default_unit=u"1", with_units=true
)
    if haskey(row, field) && !ismissing(row[field]) && haskey(row[field], :values)
        try
            val = only(row[field].values)
            if with_units
                if iszero(val)
                    val *= default_unit
                elseif haskey(row[field], :units)
                    vunit = only(get(row[field], :units, [""]))
                    val *= extract_unit(vunit, default_unit)
                else
                    val *= default_unit
                end
            end
            if verbose
                println("$(row.symbol) => $field=$val")
            end
            return val
        catch
            return missing
        end
    else
        return missing
    end
end

correct_volume_unit(v::AbstractQuantity) = uamount(v) != -1 ? v/1u"mol" : v

correct_volume_unit(v) = v

function complete_species_with_thermo_model!(species, row; verbose=false)
    Tref = row.Tst * u"K"
    Pref = row.Pst * u"Pa"
    species.Tref = Tref
    species.Pref = Pref
    values0 = [
        :Cp⁰ => extract_value(
            row, :sm_heat_capacity_p; verbose=verbose, default_unit=u"J/K/mol"
        ),
        :ΔₐH⁰ => extract_value(row, :sm_enthalpy; verbose=verbose, default_unit=u"J/mol"),
        :S⁰ =>
            extract_value(row, :sm_entropy_abs; verbose=verbose, default_unit=u"J/K/mol"),
        :ΔₐG⁰ =>
            extract_value(row, :sm_gibbs_energy; verbose=verbose, default_unit=u"J/mol"),
        :V⁰ => correct_volume_unit(extract_value(row, :sm_volume; verbose=verbose, default_unit=u"J/bar")),
    ]
    species[:thermo_params] = [values0; :T => Tref; :P => Pref]
    TPMethods = row.TPMethods
    if !ismissing(TPMethods)
        for method in TPMethods
            method_type = only(values(method.method))
            if method_type == "cp_ft_equation" && haskey(method, :m_heat_capacity_ft_coeffs)
                species[:Cp_method] = "cp_ft_equation"
                coeffs = method.m_heat_capacity_ft_coeffs
                vals = coeffs.values
                units = extract_unit.(coeffs.units)
                params = [
                    Symbol("a", subscriptnumber(i - 1)) => float(vals[i] * units[i]) for
                    i in 1:min(length(vals), length(units))
                ]
                species[:thermo_params] = [params; species[:thermo_params]]

            elseif method_type == "mv_constant"
                species[:V_method] = "mv_constant"
            end
        end
    end
    return species
end

"""
    build_species_from_database(df_substances::AbstractDataFrame, list_symbols=nothing; verbose=false) -> Vector{Species}

Build Species objects from a substance DataFrame.

# Arguments

  - `df_substances`: DataFrame containing substance data.
  - `list_symbols`: optional list of symbols to filter (default: nothing, process all).
  - `verbose`: if true, print details during processing (default: false).

# Returns

  - Vector of `Species`.
"""
function build_species_from_database(
    df_substances::AbstractDataFrame, list_symbols=nothing; verbose=false
)
    local_df_substances = if isnothing(list_symbols)
        df_substances
    else
        @view df_substances[df_substances.symbol .∈ Ref(list_symbols), :]
    end
    keylist = String[]
    species_list = Species[]
    print_title(
        "Building species"; crayon=Crayon(; foreground=:blue), style=:box, indent=""
    )
    @showprogress for row in eachrow(local_df_substances)
        if verbose
            println(row[:symbol])
        end
        species = Species(
            row.formula;
            name=row.name,
            symbol=row.symbol,
            aggregate_state=try
                eval(Meta.parse(only(values(row.aggregate_state))))
            catch
                AS_UNDEF
            end,
            class=try
                eval(Meta.parse(only(values(row.class_))))
            catch
                SC_UNDEF
            end,
        )
        complete_species_with_thermo_model!(species, row; verbose=verbose)
        key = row.symbol
        if key in keylist
            @warn("Symbol $key is used for multiple species")
        end
        push!(keylist, key)
        push!(species_list,species)
    end
    return species_list
end

function complete_reaction_with_thermo_model!(reaction, row; verbose=false)
    Tref = row.Tst * u"K"
    Pref = row.Pst * u"Pa"
    reaction.Tref = Tref
    reaction.Pref = Pref
    values0 = [
        :ΔᵣCp⁰ => extract_value(
            row, :drsm_heat_capacity_p; verbose=verbose, default_unit=u"J/K/mol"
        ),
        :ΔᵣH⁰ => extract_value(row, :drsm_enthalpy; verbose=verbose, default_unit=u"J/mol"),
        :ΔᵣS⁰ =>
            extract_value(row, :drsm_entropy_abs; verbose=verbose, default_unit=u"J/K/mol"),
        :ΔᵣG⁰ =>
            extract_value(row, :drsm_gibbs_energy; verbose=verbose, default_unit=u"J/mol"),
        :ΔᵣV⁰ => correct_volume_unit(extract_value(row, :drsm_volume; verbose=verbose, default_unit=u"J/bar")),
        :logKr => extract_value(row, :logKr; verbose=verbose, default_unit=u"1"),
    ]
    reaction[:thermo_params] = [values0; :T => Tref; :P => Pref]
    TPMethods = row.TPMethods
    if !ismissing(TPMethods)
        for method in TPMethods
            method_type = only(values(method.method))
            if method_type == "logk_fpt_function" && haskey(method, :logk_ft_coeffs)
                reaction[:logk_method] = "logk_fpt_function"
                coeffs = method.logk_ft_coeffs
                vals = coeffs.values
                units = dimension.([1, u"1/K", u"K", 1, u"K^2", u"1/K^2", u"1/√K"])
                params = [
                    Symbol("A", subscriptnumber(i - 1)) =>
                        float(Quantity(vals[i], units[i])) for
                    i in 1:min(length(vals), length(units))
                ]
                reaction[:thermo_params] = [params; reaction[:thermo_params]]
            elseif method_type == "dr_volume_constant"
                reaction[:V_method] = "dr_volume_constant"
            end
        end
    end
    return reaction
end

"""
    build_reactions_from_database(df_reactions::AbstractDataFrame, dict_species=Dict(), list_symbols=nothing; verbose=false) -> Vector{Reaction}

Build Reaction objects from a reaction DataFrame.

# Arguments

  - `df_reactions`: DataFrame containing reaction data.
  - `species_list`: vector of existing `Species` objects to use in reactions.
  - `list_symbols`: optional list of reaction symbols to filter (default: nothing, process all).
  - `verbose`: if true, print details during processing (default: false).

# Returns

  - Vector of `Reaction` objects.
"""
function build_reactions_from_database(
    df_reactions::AbstractDataFrame,
    species_list=[],
    list_symbols=nothing;
    verbose=false,
)
    local_df_reactions = if isnothing(list_symbols)
        df_reactions
    else
        @view df_reactions[df_reactions.symbol .∈ Ref(list_symbols), :]
    end
    dict_species = Dict(symbol(s) => s for s in species_list)
    keylist = String[]
    reactions_list = Reaction[]
    print_title(
        "Building reactions"; crayon=Crayon(; foreground=:red), style=:box, indent=""
    )
    function choose_species(k, rowsymbol, dict_species)
        if haskey(dict_species, k)
            return dict_species[k]
        elseif haskey(dict_species, rowsymbol)
            return dict_species[rowsymbol]
        else
            rowsymboldot = replace(rowsymbol, "_" => ".")
            if haskey(dict_species, rowsymboldot)
                return dict_species[rowsymboldot]
            else
                return find_species(k, collect(values(dict_species)))
            end
        end
    end
    @showprogress for row in eachrow(local_df_reactions)
        if verbose
            println(row[:symbol])
        end
        reaction = Reaction(
            OrderedDict(
                choose_species(last(k), row.symbol, dict_species) => last(v) for
                (k, v) in row.reactants if last(k) != "e-"
            );
            symbol=row.symbol,
        )
        complete_reaction_with_thermo_model!(reaction, row; verbose=verbose)
        key = row.symbol
        if key in keylist
            @warn("Symbol $key is used for multiple reactions")
        end
        push!(keylist, key)
        push!(reactions_list, reaction)
    end
    return reactions_list
end

"""
    get_compatible_species(df_substances::AbstractDataFrame, species_list; aggregate_states=[AS_AQUEOUS], exclude_species=[], union=false) -> DataFrame

Find species in the database compatible with a given list of species (sharing atoms).

# Arguments

  - `df_substances`: substance DataFrame.
  - `species_list`: list of target species symbols.
  - `aggregate_states`: filter for specific aggregate states (default: `[AS_AQUEOUS]`).
  - `exclude_species`: list of species symbols to exclude.
  - `union`: if true, includes the original `species_list` in the result (default: false).

# Returns

  - DataFrame of compatible substances.
"""
function get_compatible_species(
    df_substances::AbstractDataFrame,
    species_list;
    aggregate_states=[AS_AQUEOUS],
    exclude_species=[],
    union=false,
)
    df_given_species = @view df_substances[df_substances.symbol .∈ Ref(species_list), :]
    involved_atoms = union_atoms(parse_formula.(df_given_species.formula))
    mask1 = last.(only.(df_substances.aggregate_state)) .∈ Ref(string.(aggregate_states))
    mask2 = issubset.(keys.(parse_formula.(df_substances.formula)), Ref(involved_atoms))
    mask3 = .!(df_substances.symbol .∈ Ref(exclude_species))
    df_compat = @view df_substances[mask1 .&& mask2 .&& mask3, :]
    if union
        return unique(vcat(df_given_species, df_compat))
    else
        return df_compat
    end
end
