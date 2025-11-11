

abstract type Callable end  # Base type for all callable thermodynamic functions

struct ThermoFunction{P,U,F,R} <: Callable
    expr::Num
    parameters::P
    var::Num
    unit::U
    func::F
    ref::R
end

function ThermoFunction(symexpr, params_unit::NamedTuple; ref=(T=298.15u"K", P=1u"bar", t=0u"s"))
    symexpr = Num(symexpr)
    varofexpr = get_variables(symexpr)
    idxvar = findfirst(x->Symbol(x) in [:T, :t], varofexpr)
    var = Num(isnothing(idxvar) ? parse_expr_to_symbolic(:T, @__MODULE__) : varofexpr[idxvar])
    dict_vars = Dict(zip(Symbol.(varofexpr), varofexpr))
    nounitfunc = [f(var) => 1 for f in [log, log10, log2, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh, exp]]
    unit = promote_type(typeof.(values(params_unit))...) <: Quantity ? 
                Quantity(1, dimension(substitute(symexpr, [[dict_vars[k]=>v for (k,v) in pairs(params_unit)] ; nounitfunc ; var=>getfield(ref, Symbol(var))]))) : 1
    params_unit = (; (k => v for (k,v) in pairs(params_unit) if k in Symbol.(varofexpr) || isempty(varofexpr))...)
    params_nounit = (; (k=>ustrip(v) for (k,v) in pairs(params_unit))...)
    func = eval(build_function(substitute(symexpr, [dict_vars[k]=>v for (k,v) in pairs(params_nounit)]), var; expression=Val{false}))
    if iszero(symexpr) && !isempty(params_unit)
        symexpr = parse_expr_to_symbolic((keys(params_unit)[1]), @__MODULE__)
        params_unit = (; keys(params_unit)[1] => zero(values(params_unit)[1]))
    end
    return ThermoFunction(symexpr, params_unit, var, unit, func, ref)
end

function ThermoFunction(expr::Expr, params::AbstractVector{<:Number}; ref=(T=298.15u"K", P=1u"bar", t=0u"s"))
    symexpr = Num(parse_expr_to_symbolic(expr, @__MODULE__))
    varofexpr = get_variables(symexpr)
    idxvar = findfirst(x->Symbol(x) in [:T, :t], varofexpr)
    var = Num(isnothing(idxvar) ? parse_expr_to_symbolic(:T, @__MODULE__) : varofexpr[idxvar])
    deleteat!(varofexpr, findfirst(x -> isequal(x, var), varofexpr))
    sort!(varofexpr; by = x->Int(parse_symbol_num(string(x)).index))
    subszeros = Dict(k=>0 for (i,k) in enumerate(varofexpr) if length(params)<i || iszero(params[i]))
    if length(subszeros) == length(varofexpr)
        subszeros = Dict(k=>0 for (i,k) in enumerate(varofexpr) if length(params)<i || (iszero(params[i]) && i>1))
    end
    symexpr = substitute(symexpr, subszeros)
    nonzeros = (!iszero).(params)
    if !isempty(nonzeros) && !any(nonzeros)
        nonzeros[1] = true
    end
    params_unit = NamedTuple{Tuple(Symbol.(varofexpr[CartesianIndices(params)][nonzeros]))}(Tuple(params[nonzeros]))
    return ThermoFunction(symexpr, params_unit; ref=ref)
end

const thermo_function_library = Dict(
    :Cp => :(a₀ + a₁*T + a₂/T^2 + a₃/√T + a₄*T^2 + a₅*T^3 + a₆*T^4 + a₇/T^3 + a₈/T + a₉*√T + a₁₀*log(T)),
    :CpoverT => :(a₀/T + a₁ + a₂/T^3 + a₃/T^(3/2) + a₄*T + a₅*T^2 + a₆*T^3 + a₇/T^4 + a₈/T^2 + a₉/√T + a₁₀*log(T)/T),
    :logKr => :(A₀ + A₁*T + A₂/T + A₃*log(T) + A₄/T^2 + A₅*T^2 + A₆*√T),
    :∫Cp => :(((-1//2)*a₇) / (T^2) + (-a₂) / T + T*a₀ - T*a₁₀ + 2a₃*sqrt(T) + a₈*log(T) + (2//3)*(T^(3//2))*a₉ + (1//2)*(T^2)*a₁ + T*a₁₀*log(T) + (1//3)*(T^3)*a₄ + (1//4)*(T^4)*a₅ + (1//5)*(T^5)*a₆),
    :∫CpoverT => :((-2a₃) / sqrt(T) + ((-1//3)*a₇) / (T^3) + ((-1//2)*a₂) / (T^2) + (-a₈) / T + T*a₁ + a₀*log(T) + (2//1)*a₉*sqrt(T) + (1//2)*(T^2)*a₄ + (1//2)*a₁₀*(log(T)^2) + (1//3)*(T^3)*a₅ + (1//4)*(T^4)*a₆),
    :∫S => :(((1//2)*a₂) / T + ((1//6)*a₇) / (T^2) - T*a₀ + T*a₁₀ + T*a₋₁ - 4a₃*sqrt(T) - a₈*log(T) + (4//3)*(T^(3//2))*a₉ + (1//2)*(T^2)*a₁ + T*a₀*log(T) - T*a₁₀*log(T) + (1//6)*(T^3)*a₄ + (1//2)*T*a₁₀*(log(T)^2) + (1//12)*(T^4)*a₅ + (1//20)*(T^5)*a₆),
)

function ThermoFunction(sym::Symbol, params::AbstractVector{<:Number}; ref=(T=298.15u"K", P=1u"bar", t=0u"s"))
    expr = thermo_function_library[sym]
    return ThermoFunction(expr, params; ref=ref)
end

function (tf::ThermoFunction)(T)
    return tf.func(T)
end

function (tf::ThermoFunction)(T::Quantity)
    return tf.func(ustrip(T))*tf.unit
end

function (tf::ThermoFunction)()
    return tf(getfield(tf.ref, Symbol(tf.var)))
end

function get_constant_param(tf::ThermoFunction)
    indexcst = findfirst(a->iszero(expand_derivatives(ModelingToolkit.Differential(tf.var)(ModelingToolkit.Differential(parse_expr_to_symbolic(a, @__MODULE__))(tf.expr)))), keys(tf.parameters))
    if isnothing(indexcst)
        return nothing, 0
    else
        varcst = keys(tf.parameters)[indexcst]
        return varcst, expand_derivatives(ModelingToolkit.Differential(parse_expr_to_symbolic(varcst, @__MODULE__))(tf.expr))
    end
end

function +(tf::ThermoFunction, x::Number)
    varcst, α = get_constant_param(tf::ThermoFunction)
    if isnothing(varcst)
        varcst = :a₀
        if length(tf.parameters) > 0
            ps = parse_symbol_num(string(keys(tf.parameters)[1]))
            varcst = Symbol(ps.base, ps.convert_func(ps.index-1))
        end
        symvar = parse_expr_to_symbolic(varcst, @__MODULE__)
        expr = symvar + tf.expr
        return ThermoFunction(expand(expr), merge((NamedTuple{(varcst,)}(x)), tf.parameters); ref=tf.ref)
    else
        symvar = parse_expr_to_symbolic(varcst, @__MODULE__)
        expr = ModelingToolkit.expand(tf.expr + (1-α)*symvar)
        return ThermoFunction(expr, merge(tf.parameters, NamedTuple{(varcst,)}((α*getfield(tf.parameters, varcst)+x,))); ref=tf.ref)
    end
end

+(x::Number, tf::ThermoFunction) = +(tf, x)

-(tf::ThermoFunction) = ThermoFunction(-1*tf.expr, tf.parameters; ref=tf.ref)

-(tf::ThermoFunction, x::Number) = +(tf, -1*x)

-(x::Number, tf::ThermoFunction) = +(x, -1*tf)

function *(tf::ThermoFunction, x::Number)
    varcst, α = get_constant_param(tf::ThermoFunction)
    if isnothing(varcst)
        ThermoFunction(tf.expr*x, tf.parameters; ref=tf.ref)
    else
        expr = ModelingToolkit.expand(tf.expr*x + (1-α*x)*parse_expr_to_symbolic(varcst, @__MODULE__))
        return ThermoFunction(expr, merge(tf.parameters, NamedTuple{(varcst,)}((α*getfield(tf.parameters, varcst)*x,))); ref=tf.ref)
    end
end

*(x::Number, tf::ThermoFunction) = *(tf, x)

function /(tf::ThermoFunction, x::Number)
    varcst, α = get_constant_param(tf::ThermoFunction)
    if isnothing(varcst)
        ThermoFunction(tf.expr/x, tf.parameters; ref=tf.ref)
    else
        expr = ModelingToolkit.expand(tf.expr/x + (x-α)*parse_expr_to_symbolic(varcst, @__MODULE__)/x)
        return ThermoFunction(expr, merge(tf.parameters, NamedTuple{(varcst,)}((α*getfield(tf.parameters, varcst)/x,))); ref=tf.ref)
    end
end

∂(tf::ThermoFunction, var=tf.var) = ThermoFunction(expand_derivatives(ModelingToolkit.Differential(var)(tf.expr)), tf.parameters; ref=tf.ref)

∫(tf::ThermoFunction, var=tf.var) = ThermoFunction(SymbolicNumericIntegration.integrate(tf.expr, var; symbolic=true, detailed=false), tf.parameters; ref=tf.ref)

function Base.show(io::IO, tf::ThermoFunction)
    print(io, tf.expr)
    print(io, " with {")
    print(io, replace(string(tf.parameters), "("=>"", ")"=>"", " = "=>"="))
    print(io, " ; ref: ", replace(string(tf.ref), "("=>"", ")"=>"", " = "=>"="), "}")
end
