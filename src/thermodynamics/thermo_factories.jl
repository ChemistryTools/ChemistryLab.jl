const ADIM_MATH_FUNCTIONS = [:log, :log10, :log2, :log1p,
                        :exp, :expm1, :exp2, :exp10,
                        :sin, :cos, :tan, :csc, :sec, :cot,
                        :asin, :acos, :atan, :acsc, :asec, :acot,
                        :sinh, :cosh, :tanh, :csch, :sech, :coth,
                        :asinh, :acosh, :atanh, :acsch, :asech, :acoth,
                        :erf, :erfc, :gamma, :lgamma]

for f in ADIM_MATH_FUNCTIONS
    if isdefined(Base, f)
        @eval Base.$f(x::Quantity) = $f(ustrip(x))
    end
end

abstract type Callable end

"""
    ThermoFunction

Thermodynamic function with symbolic expression and compiled evaluation.
"""
struct ThermoFunction{N, R<:NamedTuple, T, D} <: Callable
    symbolic::Num
    vars::NTuple{N, Symbol}
    refs::R
    compiled::RuntimeGeneratedFunction
    unit::Quantity{T,D}
end

function extract_vars_params(expr, vars)
    params = Symbol[]
    vars_set = Set(vars)
    newvars = Symbol[]

    function scan_expr(ex)
        if ex isa Symbol
            ex ∈ vars_set ? push!(newvars, ex) : push!(params, ex)
        elseif ex isa Expr
            for arg in ex.args[2:end]
                scan_expr(arg)
            end
        end
    end

    scan_expr(expr)
    unique!(newvars)
    unique!(params)
    return newvars, params
end

"""
    ThermoFactory

Factory for creating ThermoFunctions from expressions.
"""
struct ThermoFactory
    symbolic::Num
    vars::OrderedDict{Symbol, Num}
    params::OrderedDict{Symbol, Num}
    cache::Dict{UInt64, Tuple{Num, RuntimeGeneratedFunction}}

    function ThermoFactory(symbolic::Union{SymbolicUtils.BasicSymbolic, Num}, vars::AbstractDict{Symbol, Num}, params::AbstractDict{Symbol, Num})
        new(symbolic, vars, params, Dict{UInt64, Tuple{Num, RuntimeGeneratedFunction}}())
    end
end

function ThermoFactory(expr, vars=[:T, :P, :t, :x, :y, :z]; units=nothing)
    vars, params = extract_vars_params(expr, vars)
    var_sym_dict = OrderedDict{Symbol, Num}(v => Symbolics.variable(v) for v in vars)
    param_sym_dict = OrderedDict{Symbol, Num}(p => Symbolics.variable(p) for p in params)
    to_dict(nt::NamedTuple) = Dict(pairs(nt))
    to_dict(v::AbstractVector{<:Pair}) = Dict(v)
    to_unit(s::String) = uparse(s)
    to_unit(q::Quantity) = oneunit(q)
    to_unit(x) = u"1"
    if !isnothing(units)
        dict_units = to_dict(units)
        for p in keys(var_sym_dict)
            var_sym_dict[p] = setmetadata(var_sym_dict[p], ModelingToolkitBase.VariableUnit, to_unit(dict_units[p]))
        end
        for p in keys(param_sym_dict)
            param_sym_dict[p] = setmetadata(param_sym_dict[p], ModelingToolkitBase.VariableUnit, to_unit(dict_units[p]))
        end
    end
    all_symbols = merge(var_sym_dict, param_sym_dict)
    symbolic = Symbolics.parse_expr_to_symbolic(expr, all_symbols)
    ThermoFactory(symbolic, var_sym_dict, param_sym_dict)
end

function get_unit(factory::ThermoFactory)
    return ModelingToolkitBase.get_unit(factory.symbolic)
end

function get_unit(factory::ThermoFactory, sym::Symbol)
    vp = get(factory.vars, sym, get(factory.params, sym, nothing))
    return ModelingToolkitBase.get_unit(vp)
end

"""
    (factory::ThermoFactory)(; kwargs...)

Create a ThermoFunction with caching for optimal performance.
"""
function (factory::ThermoFactory)(; kwargs...)
    param_vals = Dict{Symbol, Any}(p => get(kwargs, p, 0.0) for p in keys(factory.params))
    refs = NamedTuple([v => force_uconvert(get_unit(factory, v), kwargs[v]) for v in keys(factory.vars) if haskey(kwargs, v)])
    cache_key = hash(tuple(sort(collect(pairs(param_vals)), by=x->x.first)...))
    unit = get_unit(factory)

    simplified, compiled = get!(factory.cache, cache_key) do
        substitutions = Dict(v => safe_ustrip(get_unit(factory, p), param_vals[p]) for (p,v) in factory.params)
        substituted = Symbolics.substitute(factory.symbolic, substitutions)
        simplified = Symbolics.simplify(Symbolics.expand(substituted))
        compiled = compile_symbolic(simplified, collect(keys(factory.vars)))
        (simplified, compiled)
    end

    return ThermoFunction(simplified, Tuple(keys(factory.vars)), refs, compiled, unit)
end

"""
    compile_symbolic(symbolic_expr, var_symbols)

Compile a symbolic expression to RuntimeGeneratedFunction.
"""
function compile_symbolic(symbolic_expr, var_symbols)
    return length(var_symbols) == 1 ?
        build_function(symbolic_expr, var_symbols[1]; expression=Val(false)) :
        build_function(symbolic_expr, var_symbols...; expression=Val(false))
end

@inline function (tf::ThermoFunction{1})(; kwargs...)
    v = tf.vars[1]
    val = haskey(kwargs, v) ? kwargs[v] : tf.refs[v]
    return get(kwargs, :unit, false) ?
                tf.compiled(ustrip(val)) * tf.unit :
                tf.compiled(ustrip(val))
end

@inline function (tf::ThermoFunction{2})(; kwargs...)
    v1, v2 = tf.vars
    val1 = haskey(kwargs, v1) ? kwargs[v1] : get(tf.refs, v1, nothing)
    val2 = haskey(kwargs, v2) ? kwargs[v2] : get(tf.refs, v2, nothing)
    return get(kwargs, :unit, false) ?
                tf.compiled(ustrip(val1), ustrip(val2)) * tf.unit :
                tf.compiled(ustrip(val1), ustrip(val2))
end

@inline function (tf::ThermoFunction{3})(; kwargs...)
    v1, v2, v3 = tf.vars
    val1 = haskey(kwargs, v1) ? kwargs[v1] : get(tf.refs, v1, nothing)
    val2 = haskey(kwargs, v2) ? kwargs[v2] : get(tf.refs, v2, nothing)
    val3 = haskey(kwargs, v3) ? kwargs[v3] : get(tf.refs, v3, nothing)
    return get(kwargs, :unit, false) ?
                tf.compiled(ustrip(val1), ustrip(val2), ustrip(val3)) * tf.unit :
                tf.compiled(ustrip(val1), ustrip(val2), ustrip(val3))
end

@inline function (tf::ThermoFunction{N})(; kwargs...) where N
    if isempty(kwargs)
        var_values = ntuple(i -> ustrip(tf.refs[tf.vars[i]]), N)
        return tf.compiled(var_values...)
    else
        merged = merge(tf.refs, kwargs)
        var_values = ntuple(i -> ustrip(merged[tf.vars[i]]), N)
        return get(kwargs, :unit, false) ?
                    tf.compiled(var_values...) * tf.unit :
                    tf.compiled(var_values...)
    end
end

"""
    combine_symbolic(op, tf1::ThermoFunction, tf2::ThermoFunction)

Combine two ThermoFunctions.
"""
function combine_symbolic(op, tf1::ThermoFunction, tf2::ThermoFunction)
    all_vars = Tuple(union(tf1.vars, tf2.vars))
    refs = merge(tf1.refs, tf2.refs)

    combined = op(tf1.symbolic, tf2.symbolic)
    simplified = Symbolics.simplify(Symbolics.expand(combined))
    compiled = compile_symbolic(simplified, all_vars)

    return ThermoFunction(simplified, all_vars, refs, compiled, oneunit(op(tf1.unit, tf2.unit)))
end

"""
    combine_symbolic(op, tf::ThermoFunction, x::Number)

Combine ThermoFunction with scalar.
"""
function combine_symbolic(op, tf::ThermoFunction, x::Number)
    combined = op(tf.symbolic, ustrip(x))
    simplified = Symbolics.simplify(Symbolics.expand(combined))
    compiled = compile_symbolic(simplified, collect(tf.vars))

    return ThermoFunction(simplified, tf.vars, tf.refs, compiled, oneunit(op(tf.unit, x)))
end

"""
    combine_symbolic(op, x::Number, tf::ThermoFunction)

Combine scalar with ThermoFunction.
"""
function combine_symbolic(op, x::Number, tf::ThermoFunction)
    combined = op(ustrip(x), tf.symbolic)
    simplified = Symbolics.simplify(Symbolics.expand(combined))
    compiled = compile_symbolic(simplified, collect(tf.vars))

    return ThermoFunction(simplified, tf.vars, tf.refs, compiled, oneunit(op(x, tf.unit)))
end

"""
    apply_symbolic(op, tf::ThermoFunction)

Apply unary operation.
"""
function apply_symbolic(op, tf::ThermoFunction)
    combined = op(tf.symbolic)
    simplified = Symbolics.simplify(Symbolics.expand(combined))
    compiled = compile_symbolic(simplified, collect(tf.vars))

    return ThermoFunction(simplified, tf.vars, tf.refs, compiled, oneunit(op(tf.unit)))
end

# Binary operations
for op in (:+, :-, :*, :/, :^)
    @eval begin
        Base.$op(tf1::ThermoFunction, tf2::ThermoFunction) = combine_symbolic($op, tf1, tf2)
        Base.$op(tf::ThermoFunction, x::Number) = combine_symbolic($op, tf, x)
        Base.$op(x::Number, tf::ThermoFunction) = combine_symbolic($op, x, tf)
    end
end

Base.:-(tf::ThermoFunction) = apply_symbolic(-, tf)

for f in [ADIM_MATH_FUNCTIONS ; :sqrt; :abs]
    if isdefined(Base, f)
        @eval Base.$f(tf::ThermoFunction) = apply_symbolic($f, tf)
    end
end

function Base.show(io::IO, tf::ThermoFunction)
    print(io, tf.symbolic, " [", dimension(tf.unit), "]")
    if !isempty(tf.refs)
        print(io, " ◆ ", join(["$k=$v" for (k, v) in pairs(tf.refs)], ", "))
    end
end

function Base.show(io::IO, ::MIME"text/plain", tf::ThermoFunction)
    println(io, "ThermoFunction:")
    print(io, "  Expression: ")
    println(io, tf.symbolic, " [", dimension(tf.unit), "]")

    if !isempty(tf.refs)
        print(io, "  References: ")
        println(io, join(["$k=$v" for (k, v) in pairs(tf.refs)], ", "))
    end

    print(io, "  Variables: ")
    print(io, join(tf.vars, ", "))
end

function Base.show(io::IO, factory::ThermoFactory)
    println(io, factory.symbolic, " [", dimension(get_unit(factory)), "]")
    if !isempty(factory.params)
        print(io, " ◆ params = ", join(keys(factory.params), ", "))
    end
    if !isempty(factory.vars)
        print(io, " ◆ vars = ", join(keys(factory.vars), ", "))
    end
end

function Base.show(io::IO, ::MIME"text/plain", factory::ThermoFactory)
    println(io, "ThermoFactory:")
    print(io, "  Expression: ")
    println(io, factory.symbolic, " [", dimension(get_unit(factory)), "]")

    if !isempty(factory.params)
        print(io, "  Parameters: ")
        println(io, join(sort(collect(keys(factory.params))), ", "))
    end

    print(io, "  Variables: ")
    print(io, join(sort(collect(keys(factory.vars))), ", "))
end

# Constructors
function ThermoFunction(sym::Symbol; kwargs...)
    factory = ThermoFactory(sym, [sym]; kwargs...)
    return factory(; kwargs...)
end

function ThermoFunction(expr::Expr, vars=[:T, :P, :t, :x, :y, :z]; kwargs...)
    factory = ThermoFactory(expr, vars; kwargs...)
    return factory(; kwargs...)
end

function ThermoFunction(x::Quantity)
    x = uexpand(x)
    factory = ThermoFactory(:c; units = [:c => oneunit(x)])
    return factory(; c = x)
end

function ThermoFunction(x::Number)
    factory = ThermoFactory(:c)
    return factory(; c = x)
end

# function ThermoFunction(symexpr::Union{SymbolicUtils.BasicSymbolic, Num}; kwargs...)
#     vars = Symbol.(get_variables(symexpr))
#     refs = NamedTuple([v => kwargs[v] for v in vars if haskey(kwargs, v)])
#     compiled = compile_symbolic(symexpr, collect(vars))
#     return ThermoFunction(Num(symexpr), Tuple(vars), refs, compiled, ModelingToolkitBase.get_unit(symexpr))
# end
