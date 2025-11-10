

abstract type Callable end  # Base type for all callable thermodynamic functions

struct ThermoFunction{P,U,F,R} <: Callable
    expr::Num
    parameters::P
    var::Num
    unit::U
    func::F
    ref::R
end

function ThermoFunction(symexpr, params_unit::NamedTuple, var=nothing; ref=(T=298.15u"K", P=1u"bar", t=0u"s"))
    symexpr = Num(symexpr)
    varofexpr = Symbol.(get_variables(symexpr))
    if isnothing(var)
        vars = [:T, :t]
        idxvar = findfirst(x->x in varofexpr, vars)
        if isnothing(idxvar) idxvar = 1 end
        var = @eval (@variables $(vars[idxvar]))[1]
    end
    nounitfunc = [f(var) => 1 for f in [log, log10, log2, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh, exp]]
    unit = promote_type(typeof.(values(params_unit))...) <: Quantity ? 
                Quantity(1, dimension(substitute(symexpr, [[eval(k)=>v for (k,v) in pairs(params_unit)] ; nounitfunc ; var=>getfield(ref, Symbol(var))]))) : 1
    params_unit = (; (k => v for (k,v) in pairs(params_unit) if k in varofexpr || isempty(varofexpr))...)
    params_nounit = (; (k=>ustrip(v) for (k,v) in pairs(params_unit))...)
    func = eval(build_function(substitute(symexpr, [eval(k)=>v for (k,v) in pairs(params_nounit)]), var; expression=Val{false}))
    if iszero(symexpr) && !isempty(params_unit)
        symexpr = eval(keys(params_unit)[1])
        params_unit = (; keys(params_unit)[1] => zero(values(params_unit)[1]))
    end
    return ThermoFunction(symexpr, params_unit, var, unit, func, ref)
end

function ThermoFunction(expr::Expr, params::AbstractVector{<:Number}, var=nothing; ref=(T=298.15u"K", P=1u"bar", t=0u"s"))
    v = extract_vars(expr)
    if isnothing(var)
        vars = [:T, :t]
        idxvar = findfirst(x->x in v, vars)
        var = @eval (@variables $(vars[idxvar]))[1]
    end
    deleteat!(v, findfirst(x -> x == Symbol(var), v))
    sort!(v; by = x->Int(parse_symbol_num(string(x)).index))
    dict_vars = Dict{Symbol, Num}()
    for x in v dict_vars[x] = @eval (@parameters $x)[1] end
    subszeros = Dict(eval(k)=>0 for (i,k) in enumerate(v) if length(params)<i || iszero(params[i]))
    if length(subszeros) == length(v)
        subszeros = Dict(eval(k)=>0 for (i,k) in enumerate(v) if length(params)<i || (iszero(params[i]) && i>1))
    end
    symexpr = substitute(eval(expr), subszeros)
    nonzeros = (!iszero).(params)
    if !isempty(nonzeros) && !any(nonzeros)
        nonzeros[1] = true
    end
    params_unit = NamedTuple{Tuple(v[CartesianIndices(params)][nonzeros])}(Tuple(params[nonzeros]))
    return ThermoFunction(symexpr, params_unit, var; ref=ref)
end

const thermo_function_library = Dict(
    :Cp => (:(a₀ + a₁*T + a₂/T^2 + a₃/√T + a₄*T^2 + a₅*T^3 + a₆*T^4 + a₇/T^3 + a₈/T + a₉*√T + a₁₀*log(T)), @eval (@variables T)[1]),
    :CpoverT => (:(a₀/T + a₁ + a₂/T^3 + a₃/T^(3/2) + a₄*T + a₅*T^2 + a₆*T^3 + a₇/T^4 + a₈/T^2 + a₉/√T + a₁₀*log(T)/T), @eval (@variables T)[1]),
    :logKr => (:(A₀ + A₁*T + A₂/T + A₃*log(T) + A₄/T^2 + A₅*T^2 + A₆*√T), @eval (@variables T)[1]),
)

function ThermoFunction(sym::Symbol, params::AbstractVector{<:Number}; ref=(T=298.15u"K", P=1u"bar", t=0u"s"))
    expr, var = thermo_function_library[sym]
    return ThermoFunction(expr, params, var; ref=ref)
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
    indexcst = findfirst(a->iszero(expand_derivatives(ModelingToolkit.Differential(tf.var)(ModelingToolkit.Differential(eval(a))(tf.expr)))), keys(tf.parameters))
    if isnothing(indexcst)
        return nothing, 0
    else
        varcst = keys(tf.parameters)[indexcst]
        return varcst, expand_derivatives(ModelingToolkit.Differential(eval(varcst))(tf.expr))
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
        @eval (@parameters $(varcst))[1]
        expr = eval(varcst) + tf.expr
        return ThermoFunction(expand(expr), merge((NamedTuple{(varcst,)}(x)), tf.parameters), tf.var; ref=tf.ref)
    else
        expr = ModelingToolkit.expand(tf.expr + (1-α)*eval(varcst))
        return ThermoFunction(expr, merge(tf.parameters, NamedTuple{(varcst,)}((α*getfield(tf.parameters, varcst)+x,))), tf.var; ref=tf.ref)
    end
end

+(x::Number, tf::ThermoFunction) = +(tf, x)

-(tf::ThermoFunction) = ThermoFunction(-1*tf.expr, tf.parameters, tf.var; ref=tf.ref)

-(tf::ThermoFunction, x::Number) = +(tf, -1*x)

-(x::Number, tf::ThermoFunction) = +(x, -1*tf)

function *(tf::ThermoFunction, x::Number)
    varcst, α = get_constant_param(tf::ThermoFunction)
    if isnothing(varcst)
        ThermoFunction(tf.expr*x, tf.parameters, tf.var; ref=tf.ref)
    else
        expr = ModelingToolkit.expand(tf.expr*x + (1-α*x)*eval(varcst))
        return ThermoFunction(expr, merge(tf.parameters, NamedTuple{(varcst,)}((α*getfield(tf.parameters, varcst)*x,))), tf.var; ref=tf.ref)
    end
end

*(x::Number, tf::ThermoFunction) = *(tf, x)

function /(tf::ThermoFunction, x::Number)
    varcst, α = get_constant_param(tf::ThermoFunction)
    if isnothing(varcst)
        ThermoFunction(tf.expr/x, tf.parameters, tf.var; ref=tf.ref)
    else
        expr = ModelingToolkit.expand(tf.expr/x + (x-α)*eval(varcst)/x)
        return ThermoFunction(expr, merge(tf.parameters, NamedTuple{(varcst,)}((α*getfield(tf.parameters, varcst)/x,))), tf.var; ref=tf.ref)
    end
end

∂(tf::ThermoFunction, var=tf.var) = ThermoFunction(expand_derivatives(ModelingToolkit.Differential(var)(tf.expr)), tf.parameters, tf.var; ref=tf.ref)

∫(tf::ThermoFunction, var=tf.var) = ThermoFunction(SymbolicNumericIntegration.integrate(tf.expr, var; symbolic=true, detailed=false), tf.parameters, tf.var; ref=tf.ref)

function Base.show(io::IO, tf::ThermoFunction)
    print(io, tf.expr)
    print(io, " with {")
    print(io, replace(string(tf.parameters), "("=>"", ")"=>"", " = "=>"="))
    print(io, " ; ref: ", replace(string(tf.ref), "("=>"", ")"=>"", " = "=>"="), "}")
end
