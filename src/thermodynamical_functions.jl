

abstract type Callable end  # Base type for all callable thermodynamic functions

struct ThermoFunction{U,F,R} <: Callable
    symexpr::Num
    vars::Dict{Symbol, Num}
    unit::U
    func::F
    ref::R
end

function ThermoFunction(symexpr::Num, vars::Dict{Symbol, Num}, unit::U, ref::R) where {U,R}
    func = eval(build_function(symexpr, values(vars)...; expression=Val{false}))
    return ThermoFunction(symexpr, vars, unit, func, ref)
end

const thermo_function_library = Dict(
    :Cp => :(a₀ + a₁*T + a₂/T^2 + a₃/√T + a₄*T^2 + a₅*T^3 + a₆*T^4 + a₇/T^3 + a₈/T + a₉*√T + a₁₀*log(T)),
    :CpoverT => :(a₀/T + a₁ + a₂/T^3 + a₃/T^(3/2) + a₄*T + a₅*T^2 + a₆*T^3 + a₇/T^4 + a₈/T^2 + a₉/√T + a₁₀*log(T)/T),
    :logKr => :(A₀ + A₁*T + A₂/T + A₃*log(T) + A₄/T^2 + A₅*T^2 + A₆*√T),
)

function ThermoFunction(expr::Union{Symbol,Expr}, params = Pair[], vars = [:T, :P, :t]; ref = [:T => 298.15u"K", :P => 1u"bar", :t => 0u"s"])
    expr = get(thermo_function_library, expr, expr) 
    symexpr = Num(parse_expr_to_symbolic(expr, @__MODULE__))
    varofexpr = Num.(get_variables(symexpr))
    dictallvars = Dict(Symbol(v) => v for v ∈ varofexpr)
    dictparams = Dict(dictallvars[k] => v for (k,v) ∈ params if haskey(dictallvars, k))
    dictref = Dict(dictallvars[k] => v for (k,v) ∈ ref if haskey(dictallvars, k))
    vecvars = filter(x -> Symbol(x) ∈ vars, varofexpr)
    vecparams = filter(x -> Symbol(x) ∉ vars, varofexpr)
    veczeros = filter(x -> x ∉ keys(dictparams), vecparams)
    symexpr = substitute(symexpr, Dict(k=>0 for k ∈ veczeros))
    nounitfunc = vec([f(var) => 1 for f ∈ [log, log10, log2, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh, exp], var ∈ vecvars])
    unit = promote_type(typeof.(values(dictparams))...) <: Quantity ? 
                Quantity(1, dimension(substitute(symexpr, [collect(dictparams) ; nounitfunc ; collect(dictref)]))) : 1
    symexpr = substitute(symexpr, [k=>ustrip(v) for (k,v) in dictparams])
    dictvars = Dict{Symbol,Num}(zip(Symbol.(vecvars), vecvars))
    func = eval(build_function(symexpr, vecvars...; expression=Val{false}))
    return ThermoFunction(symexpr, dictvars, unit, func, Dict(ref))
end

function (tf::ThermoFunction)(vars...)
    return tf.func(vars...)
end

function (tf::ThermoFunction)(vars::Quantity...)
    return tf.func(ustrip.(vars)...)*tf.unit
end

function (tf::ThermoFunction)()
    vars = [tf.ref[var] for var ∈ keys(tf.vars)]
    if length(vars)>0
        return tf(vars...)
    else
        return f(0)
    end
end

function +(tf::ThermoFunction, x::Number)
    @assert dimension(tf.unit) == dimension(x)
    return ThermoFunction(tf.symexpr + ustrip(x), tf.vars, tf.unit, tf.ref)
end

+(x::Number, tf::ThermoFunction) = +(tf, x)

-(tf::ThermoFunction) = ThermoFunction(-1*tf.symexpr, tf.vars, tf.unit, tf.ref)

*(tf::ThermoFunction, x::Number) = ThermoFunction(ModelingToolkit.expand(tf.symexpr * ustrip(x)), tf.vars, tf.unit*dimension(x), tf.ref)

*(x::Number, tf::ThermoFunction) = *(tf, x)

-(tf::ThermoFunction, x) = +(tf, -1*x)

-(x::Number, tf::ThermoFunction) = +(-tf, x)

/(tf::ThermoFunction, x::Number) = ThermoFunction(ModelingToolkit.expand(tf.symexpr / ustrip(x)), tf.vars, tf.unit/dimension(x), tf.ref)

/(x::Number, tf::ThermoFunction) = ThermoFunction(ModelingToolkit.expand(ustrip(x) / tf.symexpr), tf.vars, dimension(x)/tf.unit, tf.ref)

function ^(tf::ThermoFunction, x::Number)
    @assert dimension(1) == dimension(x)
    return ThermoFunction(tf.symexpr ^ x, tf.vars, tf.unit ^ x, tf.ref)
end

function contained_dict(d1, d2)
    c1 = all(k -> haskey(d2, k) && isequal(d2[k], d1[k]), keys(d1))
    c2 = all(k -> haskey(d1, k) && isequal(d1[k], d2[k]), keys(d2))
    c1 ? d2 : c2 ? d1 : nothing
end

function +(tf1::ThermoFunction, tf2::ThermoFunction)
    vars = contained_dict(tf1.vars, tf2.vars)
    @assert !isnothing(vars) && dimension(tf1.unit) == dimension(tf2.unit)
    return ThermoFunction(tf1.symexpr + tf2.symexpr, vars, tf1.unit, merge(tf1.ref, tf2.ref))
end

function *(tf1::ThermoFunction, tf2::ThermoFunction)
    vars = contained_dict(tf1.vars, tf2.vars)
    @assert !isnothing(vars)
    return ThermoFunction(tf1.symexpr*tf2.symexpr, vars, tf1.unit*tf2.unit, merge(tf1.ref, tf2.ref))
end

function /(tf1::ThermoFunction, tf2::ThermoFunction)
    vars = contained_dict(tf1.vars, tf2.vars)
    @assert !isnothing(vars)
    return ThermoFunction(tf1.symexpr/tf2.symexpr, vars, tf1.unit/tf2.unit, merge(tf1.ref, tf2.ref))
end

Base.inv(tf::ThermoFunction) = 1/tf

∂(tf::ThermoFunction, var=collect(keys(tf.vars))[1]) = ThermoFunction(Num(ModelingToolkit.expand(expand_derivatives(ModelingToolkit.Differential(tf.vars[var])(tf.symexpr)))), tf.vars, tf.unit/dimension(tf.ref[var]), tf.ref)

∫(tf::ThermoFunction, var=collect(keys(tf.vars))[1]) = ThermoFunction(Num(ModelingToolkit.expand(SymbolicNumericIntegration.integrate(tf.symexpr, tf.vars[var]; symbolic=true, detailed=false))), tf.vars, tf.unit*dimension(tf.ref[var]), tf.ref)

function Base.show(io::IO, tf::ThermoFunction)
    print(io, tf.symexpr)
    print(io, " ♢ unit=[", dimension(tf.unit),"]")
    print(io, " ♢ ref=[")
    print(io, join(["$(k)=$(v)" for (k,v) ∈ tf.ref], ", "), "]")
end

function apply(func::Function, tf::ThermoFunction, args... ; kwargs...) 
    return ThermoFunction(Num(func(tf.symexpr, args...; kwargs...)), tf.vars, tf.unit, tf.ref)
end

for func ∈ (ModelingToolkit.expand, ModelingToolkit.simplify, sqrt, log, log10, log2, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh, exp)
    func(tf::ThermoFunction) = ThermoFunction(Num(func(tf.symexpr)), tf.vars, func(tf.unit), tf.ref)
end
