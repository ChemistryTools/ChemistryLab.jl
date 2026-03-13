"""
    THERMO_MODELS

Dictionary storing raw thermodynamic model expressions and units.
Keys are model names (symbols), values are dictionaries containing:

  - Symbolic expressions for thermodynamic functions (Cp, H, S, G).
  - Units for parameters and variables.
"""
const THERMO_MODELS = Dict(
    :cp_ft_equation => Dict(
        :Cp => :(
            a₀ +
                a₁ * T +
                a₂ / T^2 +
                a₃ / sqrt(T) +
                a₄ * T^2 +
                a₅ * T^3 +
                a₆ * T^4 +
                a₇ / T^3 +
                a₈ / T +
                a₉ * sqrt(T) +
                a₁₀ * log(T)
        ),
        :S => :(
            a₀ * log(T) +
                a₁ * T +
                -(a₂ / 2) / T^2 +
                -2 * a₃ / sqrt(T) +
                (a₄ / 2) * T^2 +
                (a₅ / 3) * T^3 +
                (a₆ / 4) * T^4 +
                -(a₇ / 3) / T^3 +
                -a₈ / T +
                2 * a₉ * sqrt(T) +
                (a₁₀ / 2) * (log(T))^2
        ),
        :H => :(
            a₀ * T +
                a₁ * T^2 / 2 +
                -a₂ / T +
                2 * a₃ * sqrt(T) +
                (a₄ / 3) * T^3 +
                (a₅ / 4) * T^4 +
                (a₆ / 5) * T^5 +
                -(a₇ / 2) / T^2 +
                a₈ * log(T) +
                (2 / 3) * a₉ * T^(3 / 2) +
                a₁₀ * T * log(T) - a₁₀ * T
        ),
        :G => :(
            -a₀ * T * log(T) +
                a₀ * T +
                -(a₁ / 2) * T^2 +
                -(a₂ / 2) / T +
                4 * a₃ * sqrt(T) +
                -(a₄ / 6) * T^3 +
                -(a₅ / 12) * T^4 +
                -(a₆ / 20) * T^5 +
                -(a₇ / 6) / T^2 +
                a₈ * log(T) +
                -(4 / 3) * a₉ * T^(3 / 2) +
                -(a₁₀ / 2) * T * (log(T))^2 +
                a₁₀ * T * log(T) - a₁₀ * T
        ),
        :units => [
            :a₀ => "J/(mol*K)",
            :a₁ => "J/(mol*K^2)",
            :a₂ => "(J*K)/mol",
            :a₃ => "J/(mol*K^0.5)",
            :a₄ => "J/(mol*K^3)",
            :a₅ => "J/(mol*K^4)",
            :a₆ => "J/(mol*K^5)",
            :a₇ => "(J*K^2)/mol",
            :a₈ => "J/mol",
            :a₉ => "J/(mol*K^1.5)",
            :a₁₀ => "J/(mol*K)",
            :T => "K",
            # :Cp => "J/(mol*K)",
            # :S => "J/(mol*K)",
            # :H => "J/mol",
            # :G => "J/mol",
        ],
    ),
    :logk_fpt_function => Dict(
        :logKr =>
            :(A₀ + A₁ * T + A₂ / T + A₃ * log(T) + A₄ / T^2 + A₅ * T^2 + A₆ * sqrt(T)),
        :units => [
            :A₀ => "1",
            :A₁ => "1/K",
            :A₂ => "K",
            :A₃ => "1",
            :A₄ => "K^2",
            :A₅ => "1/K^2",
            :A₆ => "1/√K",
            :T => "K",
            # :logKr => "1",
        ],
    ),
)

"""
    THERMO_FACTORIES

Dictionary storing compiled `ThermoFactory` objects for each model.
Used to efficiently generate `ThermoFunction` instances.
"""
const THERMO_FACTORIES = Dict{Symbol, AbstractDict}()

"""
    build_thermo_functions(model_name, params) -> OrderedDict

Build `ThermoFunction` objects for a specific model and parameters.

# Arguments

  - `model_name`: symbol identifying the thermodynamic model (e.g., `:cp_ft_equation`).
  - `params`: dictionary or pair list of parameter values.

# Returns

  - `OrderedDict` containing the constructed thermodynamic functions (`Cp⁰`, `ΔₐH⁰`, `S⁰`, `ΔₐG⁰`).
"""
function build_thermo_functions(model_name, params)
    dict_factories = THERMO_FACTORIES[model_name]
    dict_params = Dict(params)

    STref = dict_params[:S⁰]
    HTref = get(
        dict_params, :ΔfH⁰, get(dict_params, :ΔₐH⁰, get(dict_params, :ΔaH⁰, missing))
    )
    GTref = get(
        dict_params, :ΔfG⁰, get(dict_params, :ΔₐG⁰, get(dict_params, :ΔaG⁰, missing))
    )
    Tref = dict_params[:T]

    Cp⁰ = dict_factories[:Cp](; params...)

    H = dict_factories[:H](; params...)
    ΔₐH⁰ = H + (HTref - H(; T = Tref, unit = true))

    S = dict_factories[:S](; params...)
    δS⁰ = STref - S(; T = Tref, unit = true)
    S⁰ = S + δS⁰

    T = ThermoFunction(:T; units = [:T => "K"])
    if haskey(dict_factories, :G)
        G = dict_factories[:G](; params...)
        ΔₐG⁰ = (G - T * δS⁰) + (GTref - G(; T = Tref, unit = true) + Tref * δS⁰)
    else
        ΔₐG⁰ = (H - T * S⁰) + (GTref - H(; T = Tref, unit = true) + Tref * STref)
    end

    return OrderedDict(:Cp⁰ => Cp⁰, :ΔₐH⁰ => ΔₐH⁰, :S⁰ => S⁰, :ΔₐG⁰ => ΔₐG⁰)
end

"""
    build_thermo_factories(dict_expr) -> Dict

Helper function to build `ThermoFactory` objects from a model dictionary.
"""
function build_thermo_factories(dict_expr)
    return Dict(
        k => ThermoFactory(v, [:T, :P]; units = get(dict_expr, :units, nothing)) for
            (k, v) in dict_expr if k != :units
    )
end

"""
    add_thermo_model(model_name, dict_model::AbstractDict)

Add a new thermodynamic model to the registry using a dictionary of expressions.

# Arguments

  - `model_name`: unique symbol for the model.
  - `dict_model`: dictionary containing symbolic expressions and units.
"""
function add_thermo_model(model_name, dict_model::AbstractDict)
    THERMO_MODELS[model_name] = dict_model
    return THERMO_FACTORIES[model_name] = build_thermo_factories(dict_model)
end

"""
    add_thermo_model(model_name, Cpexpr::Expr, units=nothing)

Add a new thermodynamic model derived from a heat capacity (Cp) expression.
Automatically integrates Cp to find H, S, and G.

# Arguments

  - `model_name`: unique symbol for the model.
  - `Cpexpr`: symbolic expression for heat capacity as a function of T.
  - `units`: optional dictionary of units for parameters.
"""
function add_thermo_model(model_name, Cpexpr::Expr, units = nothing)
    vars, params = extract_vars_params(Cpexpr, [:T])
    var_sym_dict = Dict{Symbol, Num}(v => Symbolics.variable(v) for v in vars)
    param_sym_dict = Dict{Symbol, Num}(p => Symbolics.variable(p) for p in params)
    all_symbols = merge(var_sym_dict, param_sym_dict)
    T = var_sym_dict[:T]

    Cp = Symbolics.simplify(
        Symbolics.expand(Symbolics.parse_expr_to_symbolic(Cpexpr, all_symbols))
    )
    H = integrate(Cp, T; symbolic = true, detailed = false)

    integS = sum(terms(Cp) ./ T)
    S = integrate(integS, T; symbolic = true, detailed = false)

    G = integrate(-S, T; symbolic = true, detailed = false)

    if !isnothing(units)
        dict_units = Dict(units)
        if !haskey(dict_units, :Cp)
            push!(units, :Cp => "J/mol/K")
        end
        if !haskey(dict_units, :S)
            push!(units, :S => "J/mol/K")
        end
        if !haskey(dict_units, :H)
            push!(units, :H => "J/mol")
        end
        if !haskey(dict_units, :G)
            push!(units, :G => "J/mol")
        end
    end

    dict_model = if isnothing(units)
        Dict(:Cp => Cpexpr, :H => toexpr(H), :S => toexpr(S), :G => toexpr(G))
    else
        Dict(:Cp => Cpexpr, :H => toexpr(H), :S => toexpr(S), :G => toexpr(G), :units => units)
    end

    THERMO_MODELS[model_name] = dict_model
    return THERMO_FACTORIES[model_name] = build_thermo_factories(dict_model)
end
