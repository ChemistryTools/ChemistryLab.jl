module SymbolicNumericIntegrationExt

using ChemistryLab
import ChemistryLab:
    add_thermo_model,
    extract_vars_params,
    build_thermo_factories,
    THERMO_MODELS,
    THERMO_FACTORIES
using Symbolics: Num, variable, simplify, expand, parse_expr_to_symbolic, toexpr
using SymbolicNumericIntegration: integrate, terms

function ChemistryLab.add_thermo_model(model_name, Cpexpr::Expr, units = nothing)
    vars, params = extract_vars_params(Cpexpr, [:T])
    var_sym_dict = Dict{Symbol, Num}(v => variable(v) for v in vars)
    param_sym_dict = Dict{Symbol, Num}(p => variable(p) for p in params)
    all_symbols = merge(var_sym_dict, param_sym_dict)
    T = var_sym_dict[:T]

    Cp = simplify(expand(parse_expr_to_symbolic(Cpexpr, all_symbols)))
    H = integrate(Cp, T; symbolic = true, detailed = false)

    integS = sum(terms(Cp) ./ T)
    S = integrate(integS, T; symbolic = true, detailed = false)

    G = integrate(-S, T; symbolic = true, detailed = false)

    # Build a Dict{Symbol, Any} so we can add String default units without
    # triggering a type-mismatch push! when the caller passed Quantity values.
    units_dict = if isnothing(units)
        nothing
    else
        d = Dict{Symbol, Any}(units)
        get!(d, :Cp, "J/mol/K")
        get!(d, :S, "J/mol/K")
        get!(d, :H, "J/mol")
        get!(d, :G, "J/mol")
        collect(pairs(d))   # convert back to a vector of Pairs for ThermoFactory
    end

    dict_model = if isnothing(units_dict)
        Dict(:Cp => Cpexpr, :H => toexpr(H), :S => toexpr(S), :G => toexpr(G))
    else
        Dict(
            :Cp => Cpexpr, :H => toexpr(H), :S => toexpr(S), :G => toexpr(G),
            :units => units_dict,
        )
    end

    THERMO_MODELS[model_name] = dict_model
    return THERMO_FACTORIES[model_name] = build_thermo_factories(dict_model)
end

end # module SymbolicNumericIntegrationExt
