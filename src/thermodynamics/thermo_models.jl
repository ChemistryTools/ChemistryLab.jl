
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
        :logKr => :(A₀ + A₁ * T + A₂ / T + A₃ * log(T) + A₄ / T^2 + A₅ * T^2 + A₆ * sqrt(T)),
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

const THERMO_FACTORIES = Dict{Symbol, AbstractDict}()

function build_thermo_functions(model_name, params)
    dict_factories = THERMO_FACTORIES[model_name]
    dict_params = Dict(params)

    STref = dict_params[:S⁰]
    HTref = get(dict_params, :ΔfH⁰, get(dict_params, :ΔₐH⁰, get(dict_params, :ΔaH⁰, missing)))
    GTref = get(dict_params, :ΔfG⁰, get(dict_params, :ΔₐG⁰, get(dict_params, :ΔaG⁰, missing)))
    Tref = dict_params[:T]

    Cp⁰ = dict_factories[:Cp](; params...)

    H = dict_factories[:H](; params...)
    ΔₐH⁰ = H + (HTref - H(T = Tref; unit = true))

    S = dict_factories[:S](; params...)
    δS⁰ = STref - S(T = Tref; unit = true)
    S⁰ = S + δS⁰

    T = ThermoFunction(:T, units = [:T => "K"])
    if haskey(dict_factories, :G)
        G = dict_factories[:G](; params...)
        ΔₐG⁰ = (G - T*δS⁰) + (GTref - G(T = Tref; unit = true) + Tref*δS⁰)
    else
        ΔₐG⁰ = (H - T*S⁰) + (GTref - H(T = Tref; unit = true) + Tref*STref)
    end

    return OrderedDict(:Cp⁰ => Cp⁰, :ΔₐH⁰ => ΔₐH⁰, :S⁰ => S⁰, :ΔₐG⁰ => ΔₐG⁰)
end

function build_thermo_factories(dict_expr)
    return Dict(k => ThermoFactory(v, [:T, :P]; units=get(dict_expr, :units, nothing)) for (k, v) in dict_expr if k != :units)
end

function add_thermo_model(model_name, dict_model::AbstractDict)
    THERMO_MODELS[model_name] = dict_model
    THERMO_FACTORIES[model_name] = build_thermo_factories(dict_model)
end

function add_thermo_model(model_name, Cpexpr::Expr, units=nothing)
    vars, params = extract_vars_params(Cpexpr, [:T])
    var_sym_dict = Dict{Symbol, Num}(v => Symbolics.variable(v) for v in vars)
    param_sym_dict = Dict{Symbol, Num}(p => Symbolics.variable(p) for p in params)
    all_symbols = merge(var_sym_dict, param_sym_dict)
    T = var_sym_dict[:T]

    Cp = Symbolics.simplify(Symbolics.expand(Symbolics.parse_expr_to_symbolic(Cpexpr, all_symbols)))
    H = integrate(Cp, T; symbolic = true, detailed = false)

    integS = sum(terms(Cp)./T)
    S = integrate(integS, T; symbolic = true, detailed = false)

    G = integrate(-S, T; symbolic = true, detailed = false)

    if !isnothing(units)
        dict_units = Dict(units)
        if !haskey(dict_units, :Cp) push!(units, :Cp => "J/mol/K") end
        if !haskey(dict_units, :S) push!(units, :S => "J/mol/K") end
        if !haskey(dict_units, :H) push!(units, :H => "J/mol") end
        if !haskey(dict_units, :G) push!(units, :G => "J/mol") end
    end

    dict_model = isnothing(units) ?
                Dict(:Cp => Cpexpr, :H => toexpr(H), :S => toexpr(S), :G => toexpr(G)) :
                Dict(:Cp => Cpexpr, :H => toexpr(H), :S => toexpr(S), :G => toexpr(G), :units => units)

    THERMO_MODELS[model_name] = dict_model
    THERMO_FACTORIES[model_name] = build_thermo_factories(dict_model)
end
