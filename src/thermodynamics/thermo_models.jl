
const THERMO_MODELS = Dict(
    :cp_ft_equation => Dict(
        :Cp => :(
            a₀ +
            a₁ * T +
            a₂ / T^2 +
            a₃ / √T +
            a₄ * T^2 +
            a₅ * T^3 +
            a₆ * T^4 +
            a₇ / T^3 +
            a₈ / T +
            a₉ * √T +
            a₁₀ * log(T)
        ),
        :S => :(
            a₀ * log(T) +
            a₁ * T +
            -(a₂ / 2) / T^2 +
            -2 * a₃ / √T +
            (a₄ / 2) * T^2 +
            (a₅ / 3) * T^3 +
            (a₆ / 4) * T^4 +
            -(a₇ / 3) / T^3 +
            -a₈ / T +
            2 * a₉ * √T +
            (a₁₀ / 2) * (log(T))^2
        ),
        :H => :(
            a₀ * T +
            a₁ * T^2 / 2 +
            -a₂ / T +
            2 * a₃ * √T +
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
            4 * a₃ * √T +
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
            :Cp => "J/(mol*K)",
            :S => "J/(mol*K)",
            :H => "J/mol",
            :G => "J/mol",
        ],
    ),
    :logk_fpt_function => Dict(
        :logKr => :(A₀ + A₁ * T + A₂ / T + A₃ * log(T) + A₄ / T^2 + A₅ * T^2 + A₆ * √T),
        :units => [
            :A₀ => "1",
            :A₁ => "1/K",
            :A₂ => "K",
            :A₃ => "1",
            :A₄ => "K^2",
            :A₅ => "1/K^2",
            :A₆ => "1/√K",
            :T => "K",
            :logKr => "1",
        ],
    ),
)

const THERMO_FACTORIES = Dict{Symbol, AbstractDict}()

function check_dimensions(dict_expr, units = dict_expr[:units])
    dict_units = Dict(units)
    unit_result = Dict(k => infer_unit(v, units) for (k, v) in dict_expr if k != :units)
    for (k, v) in unit_result
        if haskey(dict_units, k)
            # println("k: ", dimension(uparse(dict_units[k])))
            # println("k: ", dimension(v))
            @assert dimension(v) === dimension(uparse(dict_units[k]))
        end
    end
    return unit_result
end

function build_thermo_functions(thermo_model_name, params)
    dict_factories = THERMO_FACTORIES[thermo_model_name]
    dict_params = Dict(params)

    STref = dict_params[:S⁰]
    HTref = get(dict_params, :ΔfH⁰, get(dict_params, :ΔₐH⁰, get(dict_params, :ΔaH⁰, missing)))
    GTref = get(dict_params, :ΔfG⁰, get(dict_params, :ΔₐG⁰, get(dict_params, :ΔaG⁰, missing)))
    Tref = dict_params[:T]

    Cp⁰ = dict_factories[:Cp](; params...)

    H = dict_factories[:H](; params...)
    ΔₐH⁰ = H + (HTref - H(T = Tref))

    S = dict_factories[:S](; params...)
    δS⁰ = STref - S(T = Tref)
    S⁰ = S + δS⁰

    T = ThermoFunction(:T)
    if haskey(dict_factories, :G)
        G = dict_factories[:G](; params...)
        ΔₐG⁰ = (G - T*δS⁰) + (GTref - G(T = Tref) + Tref*δS⁰)
    else
        ΔₐG⁰ = (H - T*S⁰) + (GTref - H(T = Tref) + Tref*STref)
    end

    return OrderedDict(:Cp⁰ => Cp⁰, :ΔₐH⁰ => ΔₐH⁰, :S⁰ => S⁰, :ΔₐG⁰ => ΔₐG⁰)
end

function build_thermo_factories(dict_expr)
    check_dimensions(dict_expr)
    return Dict(k => ThermoFactory(v, [:T, :P]) for (k, v) in dict_expr if k != :units)
end

function add_thermo_model(name, dict_model)
    THERMO_MODELS[name] = dict_model
    THERMO_FACTORIES[name] = build_thermo_factories(dict_model)
end
