using ModelingToolkit
using DynamicQuantities
using SymbolicNumericIntegration
using RuntimeGeneratedFunctions
using OrderedCollections

safe_ustrip(unit::UnionAbstractQuantity, q::UnionAbstractQuantity) = ustrip(unit, q)

safe_ustrip(::UnionAbstractQuantity, q) = ustrip(q)

safe_uconvert(qout::UnionAbstractQuantity{<:Any, <:AbstractDimensions}, q::UnionAbstractQuantity{<:Any, <:Dimensions}) = uconvert(qout, q)

safe_uconvert(::UnionAbstractQuantity{<:Any, <:AbstractDimensions}, q) = q

safe_uparse(x::AbstractString) = uparse(x)
safe_uparse(x::AbstractQuantity) = x

force_uconvert(qout::UnionAbstractQuantity{<:Any, <:AbstractDimensions}, q::UnionAbstractQuantity{<:Any, <:Dimensions}) = safe_uconvert(qout, q)
force_uconvert(qout::UnionAbstractQuantity{<:Any, <:AbstractDimensions}, q) = safe_ustrip(qout, q) * qout


const ADIM_MATH_FUNCTIONS = [
    :log, :log10, :log2, :log1p,
    :exp, :expm1, :exp2, :exp10,
    :sin, :cos, :tan, :csc, :sec, :cot,
    :asin, :acos, :atan, :acsc, :asec, :acot,
    :sinh, :cosh, :tanh, :csch, :sech, :coth,
    :asinh, :acosh, :atanh, :acsch, :asech, :acoth,
    :erf, :erfc, :gamma, :lgamma,
]

for f in ADIM_MATH_FUNCTIONS
    if isdefined(Base, f)
        @eval Base.$f(x::Quantity) = $f(ustrip(x))
    end
end

abstract type AbstractFunc end

struct SymbolicFunc{N, R <: NamedTuple, T, D} <: AbstractFunc
    symbolic::Num
    vars::NTuple{N, Symbol}
    refs::R
    compiled::RuntimeGeneratedFunction
    unit::Quantity{T, D}
end

# OPTIMIZED callable - specialized for 1, 2, 3+ variables
@inline function (tf::SymbolicFunc{1})(; kwargs...)
    # Single variable - specialized, zero allocation
    v = tf.vars[1]
    val = haskey(kwargs, v) ? kwargs[v] : tf.refs[v]
    with_unit = get(kwargs, :unit, false)
    return with_unit ? tf.compiled(ustrip(val)) * tf.unit :
        tf.compiled(ustrip(val))
end

@inline function (tf::SymbolicFunc{2})(; kwargs...)
    # Two variables - specialized
    v1, v2 = tf.vars
    val1 = haskey(kwargs, v1) ? kwargs[v1] : get(tf.refs, v1, nothing)
    val2 = haskey(kwargs, v2) ? kwargs[v2] : get(tf.refs, v2, nothing)
    with_unit = get(kwargs, :unit, false)
    return with_unit ? tf.compiled(ustrip(val1), ustrip(val2)) * tf.unit :
        tf.compiled(ustrip(val1), ustrip(val2))
end

@inline function (tf::SymbolicFunc{3})(; kwargs...)
    # Three variables - specialized
    v1, v2, v3 = tf.vars
    val1 = haskey(kwargs, v1) ? kwargs[v1] : get(tf.refs, v1, nothing)
    val2 = haskey(kwargs, v2) ? kwargs[v2] : get(tf.refs, v2, nothing)
    val3 = haskey(kwargs, v3) ? kwargs[v3] : get(tf.refs, v3, nothing)
    with_unit = get(kwargs, :unit, false)
    return with_unit ? tf.compiled(ustrip(val1), ustrip(val2), ustrip(val3)) * tf.unit :
        tf.compiled(ustrip(val1), ustrip(val2), ustrip(val3))
end

@inline function (tf::SymbolicFunc{N})(; kwargs...) where {N}
    # General case - fallback for 4+ variables
    if isempty(kwargs)
        # Fast path: no overrides, use refs directly
        var_values = ntuple(i -> ustrip(tf.refs[tf.vars[i]]), N)
        return tf.compiled(var_values...)
    else
        # Merge path
        with_unit = get(kwargs, :unit, false)
        merged = merge(tf.refs, kwargs)
        var_values = ntuple(i -> ustrip(merged[tf.vars[i]]), N)
        return with_unit ? tf.compiled(var_values...) * tf.unit : tf.compiled(var_values...)
    end
end

struct ThermoFactory
    symbolic::Num
    vars::OrderedDict{Symbol, Num}
    params::OrderedDict{Symbol, Num}
    cache::Dict{UInt64, Tuple{Num, RuntimeGeneratedFunction}}

    function ThermoFactory(symbolic::Union{SymbolicUtils.BasicSymbolic, Num}, vars::AbstractDict{Symbol, Num}, params::AbstractDict{Symbol, Num})
        return new(symbolic, vars, params, Dict{UInt64, Tuple{Num, RuntimeGeneratedFunction}}())
    end
end

function extract_vars_params(expr, vars)
    params = Symbol[]
    vars_set = Set(vars)
    newvars = Symbol[]

    function scan_expr(ex)
        return if ex isa Symbol
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

function ThermoFactory(expr, vars = [:T, :P, :t, :x, :y, :z]; units = nothing)
    vars, params = extract_vars_params(expr, vars)
    var_sym_dict = OrderedDict{Symbol, Num}(v => Symbolics.variable(v) for v in vars)
    param_sym_dict = OrderedDict{Symbol, Num}(p => Symbolics.variable(p) for p in params)
    if !isnothing(units)
        dict_units = Dict(units)
        for p in keys(var_sym_dict)
            var_sym_dict[p] = setmetadata(var_sym_dict[p], ModelingToolkitBase.VariableUnit, uparse(dict_units[p]))
        end
        for p in keys(param_sym_dict)
            param_sym_dict[p] = setmetadata(param_sym_dict[p], ModelingToolkitBase.VariableUnit, uparse(dict_units[p]))
        end
    end
    all_symbols = merge(var_sym_dict, param_sym_dict)
    symbolic = Symbolics.parse_expr_to_symbolic(expr, all_symbols)
    return ThermoFactory(symbolic, var_sym_dict, param_sym_dict)
end

function get_unit(factory::ThermoFactory)
    return ModelingToolkitBase.get_unit(factory.symbolic)
end

function get_unit(factory::ThermoFactory, sym::Symbol)
    vp = get(factory.vars, sym, get(factory.params, sym, nothing))
    return ModelingToolkitBase.get_unit(vp)
end

function compile_symbolic(symbolic_expr, var_symbols)
    return length(var_symbols) == 1 ?
        build_function(symbolic_expr, var_symbols[1]; expression = Val(false)) :
        build_function(symbolic_expr, var_symbols...; expression = Val(false))
end

function (factory::ThermoFactory)(; kwargs...)
    param_vals = Dict{Symbol, Any}(p => get(kwargs, p, 0.0) for p in keys(factory.params))
    refs = NamedTuple([v => force_uconvert(get_unit(factory, v), kwargs[v]) for v in keys(factory.vars) if haskey(kwargs, v)])
    cache_key = hash(tuple(sort(collect(pairs(param_vals)), by = x -> x.first)...))
    unit = get_unit(factory)

    simplified, compiled = get!(factory.cache, cache_key) do
        substitutions = Dict(v => param_vals[p] for (p, v) in factory.params)
        substituted = Symbolics.substitute(factory.symbolic, substitutions)
        simplified = Symbolics.simplify(Symbolics.expand(substituted))
        compiled = compile_symbolic(simplified, collect(keys(factory.vars)))

        (simplified, compiled)
    end

    return SymbolicFunc(simplified, Tuple(keys(factory.vars)), refs, compiled, unit)
end

function Base.show(io::IO, tf::SymbolicFunc)
    print(io, tf.symbolic, " [", dimension(tf.unit), "]")
    return if !isempty(tf.refs)
        print(io, " ◆ ", join(["$k=$v" for (k, v) in pairs(tf.refs)], ", "))
    end
end

function Base.show(io::IO, ::MIME"text/plain", tf::SymbolicFunc)
    println(io, "SymbolicFunc:")
    print(io, "  Expression: ")
    println(io, tf.symbolic, " [", dimension(tf.unit), "]")

    if !isempty(tf.refs)
        print(io, "  References: ")
        println(io, join(["$k=$v" for (k, v) in pairs(tf.refs)], ", "))
    end

    print(io, "  Variables: ")
    return print(io, join(tf.vars, ", "))
end

function Base.show(io::IO, factory::ThermoFactory)
    println(io, factory.symbolic, " [", dimension(get_unit(factory)), "]")
    if !isempty(factory.params)
        print(io, " ◆ params = ", join(keys(factory.params), ", "))
    end
    return if !isempty(factory.vars)
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
    return print(io, join(sort(collect(keys(factory.vars))), ", "))
end

units =
    [
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
]

dict_symvars = Dict(units)
var_sym_dict = Dict{Symbol, Num}(p => Symbolics.variable(p) for p in keys(dict_symvars))
for p in keys(dict_symvars)
    println(p)
    var_sym_dict[p] = setmetadata(var_sym_dict[p], ModelingToolkitBase.VariableUnit, uparse(dict_symvars[p]))
end

Cpexpr = :(a₀ + a₁ * T + a₂ / T^2 + a₃ / T^(1 / 2) + a₄ * T^2 + a₅ * T^3 + a₆ * T^4 + a₇ / T^3 + a₈ / T + a₉ * T^(1 / 2) + a₁₀ * log(T))
Cp = Symbolics.parse_expr_to_symbolic(Cpexpr, var_sym_dict)
ModelingToolkitBase.get_unit(Cp)
# H = integrate(Cp, var_sym_dict[:T]; symbolic = true, detailed = false)
# ModelingToolkitBase.get_unit(H)

Gexpr = :(-(a₀ * T * log(T)) + (a₀ * T) + (-(a₁ / 2) * T^2) + (-(a₂ / 2) / T) + (4 * a₃ * T^(1 / 2)) + (-(a₄ / 6) * T^3) + (-(a₅ / 12) * T^4) + (-(a₆ / 20) * T^5) + (-(a₇ / 6) / T^2) + (a₈ * log(T)) + (-(4 / 3) * a₉ * T^(3 / 2)) + (-(a₁₀ / 2) * T * (log(T))^2) + (a₁₀ * T * log(T)) + (-a₁₀ * T))
G = Symbolics.parse_expr_to_symbolic(Gexpr, var_sym_dict)
ModelingToolkitBase.get_unit(G)

factory = ThermoFactory(Cpexpr, [:T]; units = units)

ModelingToolkitBase.get_unit(factory.symbolic)
ModelingToolkitBase.get_unit(factory.params[:a₂])
ModelingToolkitBase.get_unit(factory.vars[:T])

params = [:a₀ => 210.0, :a₁ => 0.12, :a₂ => -3.07e6, :a₃ => 0.0, :T => 298.15, :Cp⁰ => 210.0, :ΔfH⁰ => -2723484.33, :S⁰ => 140, :ΔfG⁰ => -2480808.197, :V⁰ => 7.84]

tf = factory(; params...)
