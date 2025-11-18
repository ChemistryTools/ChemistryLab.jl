"""
    parse_reaction_stoich_cemdata(reaction_line::AbstractString) -> (Vector, String, String)

Parse a reaction line from CEMDATA format and extract stoichiometric information.

# Arguments

  - `reaction_line`: reaction string in CEMDATA format, optionally with a comment after '#'.

# Returns

  - `reactants`: vector of dictionaries with "symbol" and "coefficient" keys.
  - `modified_equation`: equation string with added "@" markers for aqueous species.
  - `comment`: extracted comment string (empty if none).

The function automatically adds "@" suffixes to aqueous species without explicit charges
(except for the first reactant).
"""
function parse_reaction_stoich_cemdata(reaction_line::AbstractString)
    # Split the equation and the comment
    equation_parts = split(reaction_line, '#')
    equation = strip(equation_parts[1])
    comment = length(equation_parts) > 1 ? strip(equation_parts[2]) : ""

    parts = split(equation, "=")
    if length(parts) != 2
        @warn "Equation format is incorrect: $equation"
        return [], equation, comment
    end

    modified_equation_parts = String[]
    reactants = []

    for (side_index, side) in enumerate(parts)
        tokens = split(side)
        modified_tokens = []
        for (token_index, token) in enumerate(tokens)
            if token == "+"
                push!(modified_tokens, token)
                continue
            end
            m = match(r"^([+-]?\d*\.?\d+)?(.+?)([\+\-]\d*\.?\d*)?$", token)
            if m !== nothing
                coeff_str = m.captures[1]
                base_symbol = m.captures[2]
                charge = m.captures[3]
                if isnothing(charge)
                    charge = ""
                end
                sp = base_symbol * charge
                coeff = if isnothing(coeff_str) || isempty(coeff_str)
                    1.0
                else
                    parse(Float64, coeff_str)
                end

                # Add "@" to aqueous species without charge
                if isempty(charge) && !(side_index == 1 && token_index == 1)
                    sp *= "@"
                    token = coeff_str === nothing ? sp : coeff_str * sp
                end

                push!(modified_tokens, token)

                # Add to reactants/products
                push!(
                    reactants,
                    Dict("symbol" => sp, "coefficient" => side_index == 1 ? -coeff : coeff),
                )
            else
                push!(modified_tokens, token)
            end
        end
        push!(modified_equation_parts, join(modified_tokens, " "))
    end
    modified_equation = join(modified_equation_parts, " = ")

    return reactants, modified_equation, comment
end

"""
    parse_float_array(line::AbstractString) -> Vector{Float64}

Parse a line containing space-separated floats, skipping the first token and any comments.

# Arguments

  - `line`: input string with format "keyword value1 value2 ...".

# Returns

  - Vector of successfully parsed Float64 values.

# Examples

```jldoctest
julia> parse_float_array("-analytical_expression 1.5 2.3 4.7")
3-element Vector{Float64}:
 1.5
 2.3
 4.7

julia> parse_float_array("-log_K 5.2 # comment")
1-element Vector{Float64}:
 5.2
```
"""
function parse_float_array(line)
    parts = split(line)
    if length(parts) < 2
        return Float64[]
    end
    float_parts = Float64[]
    for part in parts[2:end]
        if !startswith(part, "#")
            try
                push!(float_parts, parse(Float64, part))
            catch e
                # @warn "Could not parse '$part' as Float64, skipping."
            end
        end
    end
    return float_parts
end

"""
    parse_phases(dat_content::AbstractString) -> Dict{String,Any}

Extract phase information from PHREEQC .dat file content.

# Arguments

  - `dat_content`: full text content of a PHREEQC .dat file.

# Returns

  - Dictionary mapping phase names to their properties (equation, log_K, analytical_expression, Vm).

Parses the PHASES section and extracts reaction equations, equilibrium constants, analytical
expressions, and molar volumes for each phase.
"""
function parse_phases(dat_content)
    phases = Dict{String,Any}()
    in_phases = false
    current_phase = nothing

    for line in eachline(IOBuffer(dat_content))
        line = strip(line)
        if startswith(line, "PHASES")
            in_phases = true
            continue
        elseif in_phases && !isempty(line) && !startswith(line, "#")
            if !occursin("=", line) && !startswith(line, "-") && !isempty(line)
                parts = split(line)
                if length(parts) >= 1 && !startswith(parts[1], "-")
                    phase_name = parts[1]
                    current_phase = Dict{String,Any}("symbol" => phase_name)
                    phases[phase_name] = current_phase
                end
            elseif occursin("=", line) && current_phase !== nothing
                reactants, equation, comment = parse_reaction_stoich_cemdata(line)
                current_phase["equation"] = equation
                current_phase["reactants"] = reactants
                if !isempty(comment)
                    current_phase["comment"] = comment
                end
            elseif startswith(line, "-log_K") && current_phase !== nothing
                log_k_parts = split(line)
                if length(log_k_parts) >= 2
                    try
                        current_phase["logKr"] = Dict(
                            "values" => [parse(Float64, log_k_parts[2])], "errors" => [2]
                        )
                    catch e
                        @warn "Could not parse log_K value for phase $(current_phase["symbol"]), skipping."
                    end
                end
            elseif startswith(line, "-analytical_expression") && current_phase !== nothing
                analytical_expression = parse_float_array(line)
                # coef of log10 in .dat becomes a coef of log in .json
                if length(analytical_expression) > 3
                    analytical_expression[4] /= log(10)
                end
                current_phase["analytical_expression"] = analytical_expression
            elseif startswith(line, "-Vm") && current_phase !== nothing
                vm_parts = split(line)
                if length(vm_parts) >= 2
                    try
                        current_phase["drsm_volume"] = parse(Float64, vm_parts[2])
                    catch e
                        @warn "Could not parse Vm value for phase $(current_phase["symbol"]), skipping."
                    end
                end
            end
        end
    end

    return phases
end

"""
    merge_reactions(json_data::Dict, new_reactions::Dict) -> Dict

Merge new reactions from PHREEQC .dat into existing ThermoFun JSON data.

# Arguments

  - `json_data`: existing ThermoFun database as a dictionary.
  - `new_reactions`: dictionary of new phase reactions to add.

# Returns

  - Updated `json_data` with new reactions appended.

Only adds reactions that don't already exist (by symbol) and that have complete
required fields (logKr, analytical_expression, equation).
"""
function merge_reactions(json_data, new_reactions)
    existing_symbols = Set{String}()
    for reaction in json_data["reactions"]
        push!(existing_symbols, reaction["symbol"])
    end

    new_reactions_list = []
    for (name, phase) in new_reactions
        if !(name in existing_symbols)
            if haskey(phase, "logKr") &&
                haskey(phase, "analytical_expression") &&
                haskey(phase, "equation")
                reaction_dict = Dict{String,Any}()

                reaction_dict["symbol"] = phase["symbol"]
                reaction_dict["equation"] = phase["equation"]
                if haskey(phase, "comment")
                    reaction_dict["comment"] = phase["comment"]
                end
                reaction_dict["reactants"] = phase["reactants"]

                reaction_dict["limitsTP"] = Dict{String,Any}(
                    "range" => false,
                    "lowerP" => 0.1,
                    "lowerT" => 273.15,
                    "upperP" => 1000000,
                    "upperT" => 298.15,
                )

                reaction_dict["Tst"] = 298.15
                reaction_dict["Pst"] = 100000

                reaction_dict["TPMethods"] = [
                    Dict{String,Any}(
                        "method" => Dict("0" => "logk_fpt_function"),
                        "limitsTP" => Dict{String,Any}(
                            "lowerP" => 0,
                            "lowerT" => 273.15,
                            "upperP" => 0,
                            "upperT" => 273.15,
                        ),
                        "logk_ft_coeffs" => Dict{String,Any}(
                            "values" => vcat(
                                phase["analytical_expression"],
                                zeros(max(0, 12 - length(phase["analytical_expression"]))),
                            ),
                        ),
                    ),
                    Dict{String,Any}("method" => Dict("7" => "logk_3_term_extrap")),
                    Dict{String,Any}("method" => Dict("13" => "dr_volume_constant")),
                ]

                reaction_dict["logKr"] = phase["logKr"]

                R = 8.31446261815324
                Tst = reaction_dict["Tst"]
                logKr = phase["logKr"]["values"][1]
                dG = R * Tst * log(10) * logKr
                reaction_dict["drsm_gibbs_energy"] = Dict{String,Any}(
                    "values" => [dG], "units" => ["J/mol"]
                )

                reaction_dict["drsm_heat_capacity_p"] = Dict{String,Any}(
                    "values" => [""], "units" => ["J/(mol*K)"]
                )

                reaction_dict["drsm_enthalpy"] = Dict{String,Any}(
                    "values" => [""], "units" => ["J/mol"]
                )

                reaction_dict["drsm_entropy"] = Dict{String,Any}(
                    "values" => [""], "units" => ["J/(mol*K)"]
                )

                # The Vm field in Phreeqc .dat is the molar volume of the solid phase in cm3/mol
                # reaction_dict["drsm_volume"] = Dict{String, Any}(
                #     "values" => [phase["drsm_volume"]],
                #     "units" => ["cm3/mol"]
                # )
                reaction_dict["drsm_volume"] = Dict{String,Any}(
                    "values" => [""], "units" => ["J/bar"]
                )

                reaction_dict["datasources"] = ["Cemdata18"]

                push!(new_reactions_list, reaction_dict)
                println("New reaction added: $name")
            else
                @warn "Phase $name is missing required fields and will be skipped."
            end
        else
            println("Phase $name already exists in JSON, skipping.")
        end
    end

    append!(json_data["reactions"], new_reactions_list)
    return json_data
end

"""
    write_reaction(f::IO, reaction::Dict)

Write a single reaction dictionary to an IO stream in JSON format.

# Arguments

  - `f`: output IO stream.
  - `reaction`: reaction dictionary with all required ThermoFun fields.

Helper function used by `merge_json` to write formatted JSON output.
"""
function write_reaction(f, reaction)
    # Write the reaction as a JSON object to the file stream `f`
    write(f, "    {\n")
    write(f, "      \"symbol\": \"$(reaction["symbol"])\",\n")
    write(f, "      \"equation\": \"$(reaction["equation"])\",\n")
    if haskey(reaction, "comment")
        write(f, "      \"comment\": \"$(reaction["comment"])\",\n")
    end
    write(f, "      \"reactants\": [\n")
    for (j, reactant) in enumerate(reaction["reactants"])
        write(f, "        {\n")
        write(f, "          \"symbol\": \"$(reactant["symbol"])\",\n")
        write(f, "          \"coefficient\": $(reactant["coefficient"])\n")
        write(f, "        }")
        if j < length(reaction["reactants"])
            write(f, ",")
        end
        write(f, "\n")
    end
    write(f, "      ],\n")

    write(f, "      \"limitsTP\": {\n")
    write(f, "        \"range\": false,\n")
    write(f, "        \"lowerP\": 0.1,\n")
    write(f, "        \"lowerT\": 273.15,\n")
    write(f, "        \"upperP\": 1000000,\n")
    write(f, "        \"upperT\": 298.15\n")
    write(f, "      },\n")

    write(f, "      \"Tst\": 298.15,\n")
    write(f, "      \"Pst\": 100000,\n")

    write(f, "      \"TPMethods\": [\n")
    for (j, method) in enumerate(reaction["TPMethods"])
        write(f, "        {\n")
        if haskey(method, "method") && haskey(method["method"], "0")
            write(f, "          \"method\": {\n")
            write(f, "            \"0\": \"logk_fpt_function\"\n")
            write(f, "          },\n")
            write(f, "          \"limitsTP\": {\n")
            write(f, "            \"lowerP\": 0,\n")
            write(f, "            \"lowerT\": 273.15,\n")
            write(f, "            \"upperP\": 0,\n")
            write(f, "            \"upperT\": 273.15\n")
            write(f, "          },\n")
            write(f, "          \"logk_ft_coeffs\": {\n")
            write(f, "            \"values\": [\n")
            for (k, value) in enumerate(method["logk_ft_coeffs"]["values"])
                write(f, "              $(value)")
                if k < length(method["logk_ft_coeffs"]["values"])
                    write(f, ",")
                end
                write(f, "\n")
            end
            write(f, "            ]\n")
            write(f, "          }\n")
        elseif haskey(method, "method") && haskey(method["method"], "7")
            write(f, "          \"method\": {\n")
            write(f, "            \"7\": \"logk_3_term_extrap\"\n")
            write(f, "          }\n")
        elseif haskey(method, "method") && haskey(method["method"], "13")
            write(f, "          \"method\": {\n")
            write(f, "            \"13\": \"dr_volume_constant\"\n")
            write(f, "          }\n")
        end
        write(f, "        }")
        if j < length(reaction["TPMethods"])
            write(f, ",")
        end
        write(f, "\n")
    end
    write(f, "      ],\n")

    write(f, "      \"logKr\": {\n")
    write(f, "        \"values\": [$(reaction["logKr"]["values"][1])],\n")
    write(f, "        \"errors\": [2]\n")
    write(f, "      },\n")

    write(f, "      \"drsm_heat_capacity_p\": {\n")
    write(f, "        \"values\": [$(reaction["drsm_heat_capacity_p"]["values"][1])],\n")
    write(f, "        \"units\": [\"$(reaction["drsm_heat_capacity_p"]["values"][1])\"]\n")
    write(f, "      },\n")

    write(f, "      \"drsm_gibbs_energy\": {\n")
    write(f, "        \"values\": [$(reaction["drsm_gibbs_energy"]["values"][1])],\n")
    write(f, "        \"units\": [\"$(reaction["drsm_gibbs_energy"]["units"][1])\"]\n")
    write(f, "      },\n")

    write(f, "      \"drsm_enthalpy\": {\n")
    write(f, "        \"values\": [$(reaction["drsm_enthalpy"]["values"][1])],\n")
    write(f, "        \"units\": [\"$(reaction["drsm_enthalpy"]["units"][1])\"]\n")
    write(f, "      },\n")

    write(f, "      \"drsm_entropy\": {\n")
    write(f, "        \"values\": [$(reaction["drsm_entropy"]["values"][1])],\n")
    write(f, "        \"units\": [\"$(reaction["drsm_entropy"]["units"][1])\"]\n")
    write(f, "      },\n")

    write(f, "      \"drsm_volume\": {\n")
    write(f, "        \"values\": [$(reaction["drsm_volume"]["values"][1])],\n")
    write(f, "        \"units\": [\"$(reaction["drsm_volume"]["units"][1])\"]\n")
    write(f, "      },\n")

    write(f, "      \"datasources\": [\"Cemdata18\"]\n")

    write(f, "    }")
end

"""
    merge_json(json_path::AbstractString, dat_path::AbstractString, output_path::AbstractString)

Merge PHREEQC .dat phase data into a ThermoFun JSON database file.

# Arguments

  - `json_path`: path to input ThermoFun JSON file.
  - `dat_path`: path to PHREEQC .dat file containing phase definitions.
  - `output_path`: path for output merged JSON file.

Reads both files, extracts phases from the .dat file, merges them into the JSON
database structure, and writes the result preserving the original JSON formatting.
"""
function merge_json(json_path, dat_path, output_path)
    # Read the initial JSON file
    initial_content = read(json_path, String)

    # Parse the initial JSON file to get field order
    json_data = JSON.parsefile(json_path)

    # Preserve the initial structure
    dat_content = read(dat_path, String)
    new_reactions = parse_phases(dat_content)

    # Add new reactions
    merged_data = merge_reactions(json_data, new_reactions)

    # Write the output JSON file, preserving the initial order
    open(output_path, "w") do f
        # Find the start and end indices of the "reactions" section
        lines = split(initial_content, '\n')
        reactions_start = 0
        reactions_end = 0
        for (i, line) in enumerate(lines)
            if occursin("\"reactions\": [", line)
                reactions_start = i
            elseif reactions_start != 0 && occursin("\"elements\": [", line)
                reactions_end = i - 1
                break
            end
        end

        # Write the initial content up to the start of reactions
        for i in 1:(reactions_start - 1)
            write(f, lines[i] * "\n")
        end

        # Write the start line of reactions
        write(f, lines[reactions_start] * "\n")

        # Write all reactions (existing and new)
        for (i, reaction) in enumerate(merged_data["reactions"])
            write_reaction(f, reaction)
            if i < length(merged_data["reactions"])
                write(f, ",\n")
            else
                write(f, "\n")
            end
        end

        # Write the end of reactions and the rest of the file
        for i in reactions_end:length(lines)
            write(f, lines[i] * "\n")
        end
    end
end

"""
    extract_scal_or_vect(x) -> Any

Extract scalar from single-element vector, or return vector unchanged.

# Arguments

  - `x`: input value (potentially a vector or scalar).

# Returns

  - First element if `x` is a single-element vector, otherwise `x` unchanged.

# Examples

```jldoctest
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

"""
    extract_primary_species(file_path::AbstractString) -> DataFrame

Extract primary aqueous species from a PHREEQC database file.

# Arguments

  - `file_path`: path to PHREEQC .dat file.

# Returns

  - DataFrame with columns: species, symbol, formula, aggregate_state, atoms, charge, gamma.

Parses the SOLUTION_SPECIES section to extract master species and their properties.
The "Zz" charge placeholder is handled specially. Gamma coefficients for activity
models are extracted from "-gamma" lines.
"""
function extract_primary_species(file_path)
    lines = readlines(file_path)

    start_idx = 0
    end_idx = 0
    in_primary_section = false

    for (i, line) in enumerate(lines)
        stripped = strip(line)

        if startswith(stripped, "SOLUTION_SPECIES")
            start_idx = i + 1
            in_primary_section = true

            while start_idx <= length(lines)
                next_line = strip(lines[start_idx])
                if startswith(next_line, "# PMATCH MASTER SPECIES") ||
                    occursin("=", next_line)
                    break
                end
                start_idx += 1
            end

            if startswith(strip(lines[start_idx]), "# PMATCH MASTER SPECIES")
                start_idx += 1
            end
        end

        if in_primary_section && startswith(stripped, "# PMATCH SECONDARY MASTER SPECIES")
            end_idx = i - 1
            break
        end
    end

    if in_primary_section && end_idx == 0
        end_idx = length(lines)
    end

    species_data = []

    for i in start_idx:end_idx
        line = strip(lines[i])

        if occursin("=", line)
            parts = split(line, "=")
            current_species = strip(parts[1])

            if current_species == "e-"
                symbol = "Zz"
            elseif occursin(r"[\+\-]\d*$", current_species)
                symbol = current_species
            else
                symbol = current_species * "@"
            end

            push!(species_data, (species=current_species, symbol=symbol, gamma=Float64[]))
        end

        if startswith(line, "-gamma") && !isempty(species_data)
            parts = split(line)
            gamma_values = Float64[]
            for val in parts[2:end]
                num = tryparse(Float64, val)
                if num !== nothing
                    push!(gamma_values, num)
                end
            end

            if !isempty(gamma_values)
                last_entry = species_data[end]
                species_data[end] = (
                    species=last_entry.species, symbol=last_entry.symbol, gamma=gamma_values
                )
            end
        end
    end

    df = DataFrame(species_data)
    df.symbol = String.(df.symbol)
    df.formula .= df.symbol
    df.aggregate_state .= "AS_AQUEOUS"
    df.atoms .= parse_formula.(df.symbol)
    df.charge .= extract_charge.(df.symbol)
    df[df.symbol .== "Zz", :species] .= "Zz"
    df[df.symbol .== "Zz", :formula] .= "Zz"
    df[df.symbol .== "Zz", :charge] .= 1
    return df[sortperm(df.symbol .== "Zz"), :]
end
