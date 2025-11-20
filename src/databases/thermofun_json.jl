"""
    extract_scal_or_vect(x) -> Any

Extract scalar from single-element vector, or return vector unchanged.

# Arguments

  - `x`: input value (potentially a vector or scalar).

# Returns

  - First element if `x` is a single-element vector, otherwise `x` unchanged.

# Examples

```julia
julia> extract_scal_or_vect([42])
42

julia> extract_scal_or_vect([1, 2, 3])
3-element Vector{Int64}:
 1
 2
 3

julia> extract_scal_or_vect(missing)
missing
```
"""
function extract_scal_or_vect(x)
    if x === missing
        return x
    end
    if length(x) == 1
        return x[1]
    else
        return x
    end
end

"""
    extract_field_with_units(fields::Dict, field_name::AbstractString) -> NamedTuple

Extract value, units and errors from a ThermoFun field dictionary.

# Arguments

  - `fields`: dictionary containing thermodynamic property data.
  - `field_name`: name of the field to extract.

# Returns

  - NamedTuple with keys `values`, `units`, `errors` (each potentially `missing`).
"""
function extract_field_with_units(fields, field_name)
    if haskey(fields, field_name)
        values = get(fields[field_name], "values", missing)
        units = get(fields[field_name], "units", missing)
        errors = get(fields[field_name], "errors", missing)
        if !isa(values, Missing) && length(values) > 0
            return (
                values=extract_scal_or_vect(values),
                units=extract_scal_or_vect(units),
                errors=extract_scal_or_vect(errors),
            )
        else
            return (values=missing, units=missing, errors=missing)
        end
    else
        return (values=missing, units=missing, errors=missing)
    end
end

"""
    extract_coeffs_with_units(method::Dict, coeff_name::AbstractString) -> NamedTuple

Extract coefficient arrays with their units and names from a TPMethod.

# Arguments

  - `method`: TPMethod dictionary.
  - `coeff_name`: name of coefficient field (e.g., "logk_ft_coeffs").

# Returns

  - NamedTuple with keys `values`, `units`, `names` (each potentially `missing`).
"""
function extract_coeffs_with_units(method, coeff_name)
    if haskey(method, coeff_name)
        values = get(method[coeff_name], "values", missing)
        units = get(method[coeff_name], "units", missing)
        names = get(method[coeff_name], "names", missing)
        if !isa(values, Missing) && length(values) > 0
            return (
                values=extract_scal_or_vect(values),
                units=extract_scal_or_vect(units),
                names=extract_scal_or_vect(names),
            )
        end
    end
    return (values=missing, units=missing, names=missing)
end

"""
    extract_limitsTP(limitsTP) -> NamedTuple

Extract temperature and pressure limits from a limitsTP dictionary.

# Arguments

  - `limitsTP`: dictionary or missing.

# Returns

  - NamedTuple with keys `lowerT`, `upperT`, `lowerP`, `upperP`, `range`.
"""
function extract_limitsTP(limitsTP)
    if limitsTP === missing
        return (
            lowerT=missing, upperT=missing, lowerP=missing, upperP=missing, range=missing
        )
    else
        return (
            lowerT=get(limitsTP, "lowerT", missing),
            upperT=get(limitsTP, "upperT", missing),
            lowerP=get(limitsTP, "lowerP", missing),
            upperP=get(limitsTP, "upperP", missing),
            range=get(limitsTP, "range", missing),
        )
    end
end

"""
    parse_TPMethod(method::Dict) -> Dict

Parse a single TPMethod entry from ThermoFun format.

# Arguments

  - `method`: raw TPMethod dictionary.

# Returns

  - Parsed dictionary with method name, limitsTP, and coefficient arrays.
"""
function parse_TPMethod(method)
    method_dict = Dict()
    method_dict["method"] = only(values(get(method, "method", Dict("" => missing))))
    if haskey(method, "limitsTP")
        method_dict["limitsTP"] = extract_limitsTP(method["limitsTP"])
    end
    if haskey(method, "logk_ft_coeffs")
        method_dict["logk_ft_coeffs"] = extract_coeffs_with_units(method, "logk_ft_coeffs")
    end
    if haskey(method, "m_heat_capacity_ft_coeffs")
        method_dict["m_heat_capacity_ft_coeffs"] = extract_coeffs_with_units(
            method, "m_heat_capacity_ft_coeffs"
        )
    end
    return method_dict
end

"""
    get_value(row, field::Symbol; debug=false, crayon=Crayon(), with_units=true, default_unit=u"1") -> Any

Extract a thermodynamic property value from a DataFrame row with optional unit attachment.

# Arguments

  - `row`: DataFrame row containing substance or reaction data.
  - `field`: property field name to extract.
  - `debug`: if > 1, print debug information for zero or missing values.
  - `crayon`: Crayon for colored debug output.
  - `with_units`: if true, attach units using DynamicQuantities.
  - `default_unit`: fallback unit if parsing fails.

# Returns

  - Property value with or without units, or `missing` if field absent.
"""
function get_value(
    row, field::Symbol; debug=false, crayon=Crayon(), with_units=true, default_unit=u"1"
)
    if haskey(row, field)
        val = row[field].values
        vunit = row[field].units
        if debug > 1 && (iszero(val) || (with_units && ismissing(vunit)))
            println(crayon("$(row.symbol) => $field=$val $vunit"))
            print(crayon"reset")
        end
        if with_units
            val = val * (
                try
                    uparse(vunit)
                catch
                    default_unit
                end
            )
        end
        return val
    else
        return missing
    end
end

"""
    get_Cp_coef(row; debug=false, crayon=Crayon(), with_units=true) -> Union{Vector,Missing}

Extract heat capacity coefficients from a substance row.

# Arguments

  - `row`: DataFrame row with substance data.
  - `debug`: if > 1, print debug info for non-zero higher-order terms.
  - `crayon`: Crayon for colored debug output.
  - `with_units`: if true, return coefficients as `Symbol => Quantity` pairs.

# Returns

  - Vector of coefficient pairs `[:a₀ => value, :a₁ => value, ...]` or `missing`.

Looks for "cp_ft_equation" TPMethod; falls back to constant Cp⁰ if unavailable.
"""
function get_Cp_coef(row; debug=false, crayon=Crayon(), with_units=true)
    TPMethods = row.TPMethods
    idx = findfirst(d -> haskey(d, "method") && d["method"] == "cp_ft_equation", TPMethods)
    if !isnothing(idx)
        d = TPMethods[idx]
        tuple_coefs = d["m_heat_capacity_ft_coeffs"]
        values = tuple_coefs.values
        if debug > 1 && !iszero(max(abs.(values[5:end])...))
            println(crayon("$(row.symbol) => Cp=$values"))
        end
        if with_units
            units = tuple_coefs.units
            values = [
                Symbol("a", subscriptnumber(i - 1)) => float(values[i] * uparse(units[i]))
                for i in 1:min(length(values), length(units))
            ]
        end
        return values
    else
        Cp = get_value(
            row,
            :Cp⁰;
            debug=debug,
            crayon=crayon,
            with_units=with_units,
            default_unit=u"J/(mol*K)",
        )
        if !ismissing(Cp)
            return [:a₀ => Cp]
        else
            return missing
        end
    end
end

"""
    complete_species_database!(df_substances::DataFrame; with_units=true, debug=false, all_properties=false) -> DataFrame

Complete a substances DataFrame by constructing Species objects with thermodynamic functions.

# Arguments

  - `df_substances`: DataFrame with ThermoFun substance data.
  - `with_units`: attach DynamicQuantities units to properties.
  - `debug`: enable debug output.
  - `all_properties`: if true, compute temperature-dependent functions (Cp⁰, ΔfH⁰, S⁰, ΔfG⁰).

# Returns

  - Modified DataFrame with added `:species` and `:cemspecies` columns.

Constructs `Species` objects from formula and properties, builds `ThermoFunction`
instances for Cp⁰(T), ΔfH⁰(T), S⁰(T), ΔfG⁰(T) when `all_properties=true`.
"""
function complete_species_database!(
    df_substances::DataFrame; with_units=true, debug=false, all_properties=false
)
    colspecies = [
        try
            Species(
                s.formula;
                name=s.name,
                symbol=s.symbol,
                aggregate_state=try
                    eval(Meta.parse(s.aggregate_state))
                catch
                    AS_UNDEF
                end,
                class=try
                    eval(Meta.parse(s.class))
                catch
                    SC_UNDEF
                end,
            )
        catch
            missing
        end for s in eachrow(df_substances)
    ]

    insertcols!(df_substances, 1, :species => colspecies)

    print_title(
        "Species property completion";
        crayon=Crayon(; foreground=:blue),
        style=:box,
        indent="",
    )

    @showprogress for row in eachrow(df_substances)
        s = row.species
        if debug
            println(s)
        end
        Tref = with_units ? row.Tst * u"K" : row.Tst
        Pref = with_units ? row.Pst * u"Pa" : row.Pst
        s.Tref = Tref
        s.Pref = Pref
        if !with_units
            s.M = ustrip(s.M)
        end

        if all_properties
            coeffa = get_Cp_coef(
                row; debug=debug, crayon=crayon"green", with_units=with_units
            )

            if !ismissing(coeffa)
                s.Cp⁰ = ThermoFunction(:Cp, coeffa; ref=[:T => Tref, :P => Pref])

                ΔfH⁰ = get_value(
                    row,
                    :ΔfH⁰;
                    debug=debug,
                    crayon=crayon"red",
                    with_units=with_units,
                    default_unit=u"J/mol",
                )
                ∫Cp = ∫(s.Cp⁰)
                s.ΔfH⁰ = (ΔfH⁰ - ∫Cp(Tref)) + ∫Cp

                S0 = get_value(
                    row,
                    :S⁰;
                    debug=debug,
                    crayon=crayon"red",
                    with_units=with_units,
                    default_unit=u"J/(mol*K)",
                )
                CpoverT = ThermoFunction(:CpoverT, coeffa; ref=[:T => Tref, :P => Pref])
                ∫CpoverT = ∫(CpoverT)
                s.S⁰ = (S0 - ∫CpoverT(Tref)) + ∫CpoverT

                ΔfG⁰ = get_value(
                    row,
                    :ΔfG⁰;
                    debug=debug,
                    crayon=crayon"blue",
                    with_units=with_units,
                    default_unit=u"J/mol",
                )
                ∫S = ∫(s.S⁰)
                s.ΔfG⁰ = (ΔfG⁰ + ∫S(Tref)) - ∫S

                Vm = get_value(
                    row, :Vm; crayon=crayon"blue", with_units=true, default_unit=u"J/bar"
                )
                s.Vm = with_units ? Vm / u"mol" : ustrip(Vm)
            end
        end
    end

    colcemspecies = [
        try
            CemSpecies(s)
        catch
            missing
        end for s in colspecies
    ]
    insertcols!(df_substances, 2, :cemspecies => colcemspecies)

    return df_substances
end

"""
    get_logKr_coef(row; debug=false, crayon=Crayon(), with_units=true) -> Union{Vector,Missing}

Extract log K equilibrium constant coefficients from a reaction row.

# Arguments

  - `row`: DataFrame row with reaction data.
  - `debug`: if > 1, print debug info for non-zero higher-order terms.
  - `crayon`: Crayon for colored debug output.
  - `with_units`: if true, return coefficients as `Symbol => Quantity` pairs.

# Returns

  - Vector of coefficient pairs `[:A₀ => value, :A₁ => value, ...]` or `missing`.

Looks for "logk_fpt_function" TPMethod; falls back to constant logKr if unavailable.
"""
function get_logKr_coef(row; debug=false, crayon=Crayon(), with_units=true)
    TPMethods = row.TPMethods
    idx = findfirst(
        d -> haskey(d, "method") && d["method"] == "logk_fpt_function", TPMethods
    )
    if !isnothing(idx)
        d = TPMethods[idx]
        if haskey(d, "logk_ft_coeffs")
            tuple_coefs = d["logk_ft_coeffs"]
            values = tuple_coefs.values
            if debug > 1 && !iszero(max(abs.(values[8:end])...))
                println(crayon("$(row.symbol) => logKr=$values"))
            end
            if with_units
                units = dimension.([1, u"1/K", u"K", 1, u"K^2", u"1/K^2", u"1/√K"])
                # values = [values[i]*units[i] for i=1:min(length(values), length(units))]
                values = [
                    Symbol("A", subscriptnumber(i - 1)) =>
                        float(Quantity(values[i], units[i])) for
                    i in 1:min(length(values), length(units))
                ]
            end
            return values
        else
            return missing
        end
    else
        logKr = get_value(row, :logKr; debug=debug, crayon=crayon, with_units=false)
        if !ismissing(logKr)
            return [:A₀ => logKr]
        else
            return missing
        end
    end
end

"""
    complete_reaction_database!(df_reactions::DataFrame, df_substances=nothing; with_units=true, debug=false, all_properties=false) -> DataFrame

Complete a reactions DataFrame by constructing Reaction objects with thermodynamic functions.

# Arguments

  - `df_reactions`: DataFrame with ThermoFun reaction data.
  - `df_substances`: optional DataFrame with substance data for matching species.
  - `with_units`: attach DynamicQuantities units to properties.
  - `debug`: enable debug output.
  - `all_properties`: if true, populate temperature-dependent logKr and reaction properties.

# Returns

  - Modified DataFrame with added `:reaction` column.

Constructs `Reaction` objects from reactant dictionaries, builds `ThermoFunction`
instances for logKr(T) when `all_properties=true`.
"""
function complete_reaction_database!(
    df_reactions::DataFrame,
    df_substances=nothing;
    with_units=true,
    debug=false,
    all_properties=false,
)
    print_title(
        "Reaction property completion";
        crayon=Crayon(; foreground=:red),
        style=:box,
        indent="",
    )

    function populate_reaction(r, row)
        if debug
            println(r)
        end
        Tref = with_units ? row.Tst * u"K" : row.Tst
        Pref = with_units ? row.Pst * u"Pa" : row.Pst
        r.Tref = Tref
        r.Pref = Pref

        if all_properties
            coefflogKr = get_logKr_coef(
                row; debug=debug, crayon=crayon"green", with_units=with_units
            )

            if !ismissing(coefflogKr)
                r.logKr = ThermoFunction(:logKr, coefflogKr; ref=[:T => Tref, :P => Pref])
            end

            r.ΔrCp⁰_Tref = get_value(
                row,
                :ΔrCp⁰;
                debug=debug,
                crayon=crayon"red",
                with_units=with_units,
                default_unit=u"J/(mol*K)",
            )
            r.ΔrH⁰_Tref = get_value(
                row,
                :ΔrH⁰;
                debug=debug,
                crayon=crayon"red",
                with_units=with_units,
                default_unit=u"J/mol",
            )
            r.ΔrG⁰_Tref = get_value(
                row,
                :ΔrG⁰;
                debug=debug,
                crayon=crayon"red",
                with_units=with_units,
                default_unit=u"J/mol",
            )
            r.ΔrS⁰_Tref = get_value(
                row,
                :ΔrS⁰;
                debug=debug,
                crayon=crayon"red",
                with_units=with_units,
                default_unit=u"J/(mol*K)",
            )
            ΔVm = get_value(
                row,
                :ΔrV;
                debug=debug,
                crayon=crayon"red",
                with_units=true,
                default_unit=u"J/bar",
            )
            r.ΔrV_Tref = with_units ? ΔVm / u"mol" : ustrip(ΔVm)
            r.logKr_Tref = get_value(
                row,
                :logKr;
                debug=debug,
                crayon=crayon"red",
                with_units=with_units,
                default_unit=u"1",
            )
        end
        return r
    end

    function species_or_row_symbol(k, reaction_symbol)
        if k == reaction_symbol
            return k
        end
        if !isnothing(df_substances) &&
            size(
            filter(
                x ->
                    x.symbol == reaction_symbol &&
                    (x.formula == k || occursin(k, x.symbol)),
                df_substances,
            ),
            1,
        ) == 1
            return reaction_symbol
        else
            return k
        end
    end

    species_list = if isnothing(df_substances)
        nothing
    elseif "species" in names(df_substances)
        df_substances.species
    else
        nothing
    end

    reactions = @showprogress [
        populate_reaction(
            Reaction(
                Dict(
                    find_species(species_or_row_symbol(k, row.symbol), species_list) => v
                    for (k, v) in row.reactants if k != "e-"
                ),
            ),
            row,
        ) for row in eachrow(df_reactions)
    ]
    insertcols!(df_reactions, :reaction => reactions)

    return df_reactions
end

"""
    read_thermofun_substances(substances::JSON3.Array; with_units=true, add_species=true, all_properties=true, debug=false) -> DataFrame
    read_thermofun_substances(filename::AbstractString; with_units=true, add_species=true, all_properties=true, debug=false) -> DataFrame

Parse ThermoFun substances data into a DataFrame.

# Arguments

  - `substances`: JSON3 array of substance records, or
  - `filename`: path to ThermoFun JSON file.
  - `with_units`: attach DynamicQuantities units.
  - `add_species`: if true, construct Species objects.
  - `all_properties`: if true, build temperature-dependent thermodynamic functions.
  - `debug`: enable debug output.

# Returns

  - DataFrame with substance properties and optionally Species objects.
"""
function read_thermofun_substances(
    substances::JSON3.Array;
    with_units=true,
    add_species=true,
    all_properties=true,
    debug=false,
)
    df_substances = DataFrame(;
        name=[get(s, "name", missing) for s in substances],
        symbol=[get(s, "symbol", missing) for s in substances],
        formula=[get(s, "formula", missing) for s in substances],
        charge=[get(s, "formula_charge", missing) for s in substances],
        aggregate_state=[
            only(values(get(s, "aggregate_state", Dict("" => missing)))) for s in substances
        ],
        class=[only(values(get(s, "class_", Dict("" => missing)))) for s in substances],
        Tst=[get(s, "Tst", missing) for s in substances],
        Pst=[get(s, "Pst", missing) for s in substances],
        TPMethods=[
            [parse_TPMethod(m) for m in get(s, "TPMethods", [])] for s in substances
        ],
        Cp⁰=[extract_field_with_units(s, "sm_heat_capacity_p") for s in substances],
        ΔfG⁰=[extract_field_with_units(s, "sm_gibbs_energy") for s in substances],
        ΔfH⁰=[extract_field_with_units(s, "sm_enthalpy") for s in substances],
        S⁰=[extract_field_with_units(s, "sm_entropy_abs") for s in substances],
        Vm=[extract_field_with_units(s, "sm_volume") for s in substances],
        datasources=[get(s, "datasources", missing) for s in substances],
    )
    if add_species
        complete_species_database!(
            df_substances; with_units=with_units, all_properties=all_properties, debug=debug
        )
    end

    return df_substances
end

function read_thermofun_substances(
    filename::AbstractString;
    with_units=true,
    add_species=true,
    all_properties=true,
    debug=false,
)
    read_thermofun_substances(
        open(JSON3.read, filename).substances;
        with_units=with_units,
        add_species=add_species,
        all_properties=all_properties,
        debug=debug,
    )
end

"""
    read_thermofun_reactions(reactions::JSON3.Array, df_substances=nothing; with_units=true, add_reactions=true, all_properties=true, debug=false) -> DataFrame
    read_thermofun_reactions(filename::AbstractString, df_substances=nothing; with_units=true, add_reactions=true, all_properties=true, debug=false) -> DataFrame

Parse ThermoFun reactions data into a DataFrame.

# Arguments

  - `reactions`: JSON3 array of reaction records, or
  - `filename`: path to ThermoFun JSON file.
  - `df_substances`: optional DataFrame with substances for species matching.
  - `with_units`: attach DynamicQuantities units.
  - `add_reactions`: if true, construct Reaction objects.
  - `all_properties`: if true, build temperature-dependent logKr functions.
  - `debug`: enable debug output.

# Returns

  - DataFrame with reaction properties and optionally Reaction objects.
"""
function read_thermofun_reactions(
    reactions::JSON3.Array,
    df_substances=nothing;
    with_units=true,
    add_reactions=true,
    all_properties=true,
    debug=false,
)
    df_reactions = DataFrame(;
        symbol=[get(r, "symbol", missing) for r in reactions],
        equation=[get(r, "equation", missing) for r in reactions],
        reactants=[
            Dict(r["symbol"] => r["coefficient"] for r in get(r, "reactants", [])) for
            r in reactions
        ],
        limitsTP=[extract_limitsTP(get(r, "limitsTP", missing)) for r in reactions],
        Tst=[get(r, "Tst", missing) for r in reactions],
        Pst=[get(r, "Pst", missing) for r in reactions],
        TPMethods=[[parse_TPMethod(m) for m in get(r, "TPMethods", [])] for r in reactions],
        logKr=[extract_field_with_units(r, "logKr") for r in reactions],
        ΔrCp⁰=[extract_field_with_units(r, "drsm_heat_capacity_p") for r in reactions],
        ΔrG⁰=[extract_field_with_units(r, "drsm_gibbs_energy") for r in reactions],
        ΔrH⁰=[extract_field_with_units(r, "drsm_enthalpy") for r in reactions],
        ΔrS⁰=[extract_field_with_units(r, "drsm_entropy") for r in reactions],
        ΔrV=[extract_field_with_units(r, "drsm_volume") for r in reactions],
        datasources=[get(r, "datasources", missing) for r in reactions],
    )
    if add_reactions
        complete_reaction_database!(
            df_reactions,
            df_substances;
            with_units=with_units,
            debug=debug,
            all_properties=all_properties,
        )
    end

    return df_reactions
end

function read_thermofun_reactions(
    filename::AbstractString,
    df_substances=nothing;
    with_units=true,
    add_reactions=true,
    all_properties=true,
    debug=false,
)
    read_thermofun_reactions(
        open(JSON3.read, filename).reactions,
        df_substances;
        with_units=with_units,
        add_reactions=add_reactions,
        all_properties=all_properties,
        debug=debug,
    )
end

"""
    read_thermofun_elements(elements::JSON3.Array) -> DataFrame
    read_thermofun_elements(filename::AbstractString) -> DataFrame

Parse ThermoFun elements data into a DataFrame.

# Arguments

  - `elements`: JSON3 array of element records, or
  - `filename`: path to ThermoFun JSON file.

# Returns

  - DataFrame with element properties (symbol, class, entropy, atomic mass, datasources).
"""
function read_thermofun_elements(elements::JSON3.Array)
    df_elements = DataFrame(;
        symbol=[get(e, "symbol", missing) for e in elements],
        class=[only(values(get(e, "class_", Dict("" => missing)))) for e in elements],
        S=[extract_field_with_units(e, "entropy") for e in elements],
        atomic_mass=[extract_field_with_units(e, "atomic_mass") for e in elements],
        datasources=[get(e, "datasources", missing) for e in elements],
    )

    return df_elements
end

function read_thermofun_elements(filename::AbstractString)
    read_thermofun_elements(open(JSON3.read, filename).elements)
end

"""
    read_thermofun(filename::AbstractString; with_units=true, add_species=true, add_reactions=true, all_properties=true, debug=false) -> Tuple{DataFrame, DataFrame, DataFrame}

Read a complete ThermoFun database file and parse all sections.

# Arguments

  - `filename`: path to ThermoFun JSON database file.
  - `with_units`: attach DynamicQuantities units.
  - `add_species`: construct Species objects.
  - `add_reactions`: construct Reaction objects.
  - `all_properties`: build temperature-dependent thermodynamic functions.
  - `debug`: enable debug output.

# Returns

  - Tuple of three DataFrames: `(df_elements, df_substances, df_reactions)`.
"""
function read_thermofun(
    filename;
    with_units=true,
    add_species=true,
    add_reactions=true,
    all_properties=true,
    debug=false,
)
    data = open(JSON3.read, filename)
    df_elements = read_thermofun_elements(data.elements)
    df_substances = read_thermofun_substances(
        data.substances;
        with_units=with_units,
        add_species=add_species,
        all_properties=all_properties,
        debug=debug,
    )
    df_reactions = read_thermofun_reactions(
        data.reactions,
        df_substances;
        with_units=with_units,
        add_reactions=add_reactions,
        all_properties=all_properties,
        debug=debug,
    )
    return df_elements, df_substances, df_reactions
end
