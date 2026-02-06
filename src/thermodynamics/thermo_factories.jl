using ModelingToolkit
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

"""
    calculate_molar_mass(atoms::AbstractDict{Symbol,T}) where {T<:Number} -> Quantity

Calculate the molar mass from an atomic composition dictionary.

# Arguments

  - `atoms`: dictionary mapping element symbols to stoichiometric coefficients.

# Returns

  - Molar mass as a Quantity in g/mol units.

# Examples

```jldoctest
julia> calculate_molar_mass(OrderedDict(:H => 2, :O => 1))
0.0180149999937744 kg mol⁻¹
```
"""
function calculate_molar_mass(atoms::AbstractDict{Symbol,T}) where {T<:Number}
    # return sum(cnt * ustrip(elements[element].atomic_mass) for (element, cnt) in atoms if haskey(elements, element); init=0) * u"g/mol"
    # return uconvert(u"g/mol", sum(cnt * elements[element].atomic_mass for (element, cnt) in atoms if haskey(elements, element); init=0u) * AvogadroConstant)
    molar_masses = [
        cnt * convert(DynamicQuantities.Quantity, elements[element].atomic_mass) for
        (element, cnt) in atoms if haskey(elements, element)
    ]
    return length(molar_masses) > 0 ? sum(molar_masses) * Constants.N_A : 0u"g/mol"
end

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

"""
    infer_unit(expr::Expr, params::Vector)

Infer the output unit of a symbolic expression while checking dimensional consistency.

Transcendental functions (log, exp, sin, cos, etc.) accept dimensional arguments
and always return dimensionless quantities.

```jldoctest
julia> expr = :(α + β * T + γ * log(T))
:(α + β * T + γ * log(T))

julia> params = [:α => "J/K/mol", :β => "J/mol/K^2", :γ => "J/K/mol"]
3-element Vector{Pair{Symbol, String}}:
 :α => "J/K/mol"
 :β => "J/mol/K^2"
 :γ => "J/K/mol"

julia> ref = [:T => "K"]
1-element Vector{Pair{Symbol, String}}:
 :T => "K"

julia> infer_unit(expr, params, ref)
1.0 m² kg s⁻² K⁻¹ mol⁻¹

julia> infer_unit(expr, params, ref; unit=us"J/mol/K")
1.0 K⁻¹ mol⁻¹ J
```
"""
function infer_unit(expr::Expr, params, ref=[]; unit=nothing)
    param_dict = Dict(kv.first => safe_uparse(kv.second) for kv in params)
    for kv in ref
        param_dict[kv.first] = safe_uparse(kv.second)
    end
    result = _infer_unit(expr, param_dict)
    return isnothing(unit) ? result : uconvert(unit, result)
end

function _infer_unit(expr, param_dict::Dict)
    if expr isa Symbol
        # Symbol: retrieve its unit
        if haskey(param_dict, expr)
            return param_dict[expr]
        else
            throw(ArgumentError("Symbol $expr not defined in parameters"))
        end

    elseif expr isa Number
        # Number: dimensionless (create a dimensionless quantity)
        return Quantity(float(expr), dimension=NoDims)

    elseif expr isa Expr
        if expr.head == :call
            func = expr.args[1]
            args = expr.args[2:end]

            # Transcendental functions: arguments can have dimensions, output is dimensionless
            if func in ADIM_MATH_FUNCTIONS
                # Compute units of arguments (to verify they are well-defined)
                # but don't check their dimension
                for arg in args
                    _infer_unit(arg, param_dict)
                end
                return Quantity(1.0, dimension=NoDims)  # Return dimensionless

            # Arithmetic operations
            elseif func == :+
                units = [_infer_unit(arg, param_dict) for arg in args]
                # Check that all units are compatible
                ref_unit = units[1]
                for u in units[2:end]
                    if dimension(u) != dimension(ref_unit)
                        throw(DimensionError("Unit incompatibility in addition: $(dimension(ref_unit)) vs $(dimension(u))"))
                    end
                end
                return ref_unit

            elseif func == :-
                if length(args) == 1
                    # Unary negation
                    return _infer_unit(args[1], param_dict)
                else
                    # Subtraction
                    units = [_infer_unit(arg, param_dict) for arg in args]
                    ref_unit = units[1]
                    for u in units[2:end]
                        if dimension(u) != dimension(ref_unit)
                            throw(DimensionError("Unit incompatibility in subtraction: $(dimension(ref_unit)) vs $(dimension(u))"))
                        end
                    end
                    return ref_unit
                end

            elseif func == :*
                units = [_infer_unit(arg, param_dict) for arg in args]
                result = units[1]
                for u in units[2:end]
                    result = result * u
                end
                return result

            elseif func == :/
                unit_num = _infer_unit(args[1], param_dict)
                unit_den = _infer_unit(args[2], param_dict)
                return unit_num / unit_den

            elseif func == :^
                base_unit = _infer_unit(args[1], param_dict)
                unit_exponent = _infer_unit(args[2], param_dict)
                return base_unit ^ unit_exponent
                # exponent = args[2]

                # # Exponent must be a dimensionless number
                # if exponent isa Number
                #     return base_unit ^ exponent
                # else
                #     throw(ArgumentError("Exponent in power must be a number"))
                # end

            else
                base_unit = _infer_unit(args[1], param_dict)
                return eval(func)(base_unit)
                # throw(ArgumentError("Unsupported function: $func"))
            end
        else
            throw(ArgumentError("Unsupported expression: $(expr.head)"))
        end
    else
        throw(ArgumentError("Unsupported expression type: $(typeof(expr))"))
    end
end

abstract type Callable end

"""
    ThermoFunction

Thermodynamic function with symbolic expression and compiled evaluation.
"""
struct ThermoFunction{N, R<:NamedTuple} <: Callable
    symbolic::Num
    vars::NTuple{N, Symbol}
    refs::R
    compiled::RuntimeGeneratedFunction
end

"""
    ThermoFactory

Factory for creating ThermoFunctions from expressions.
"""
struct ThermoFactory{E}
    expr::E
    vars::Vector{Symbol}
    params::Vector{Symbol}
    cache::Dict{UInt64, Tuple{Num, RuntimeGeneratedFunction}}

    function ThermoFactory{E}(expr::E, vars::Vector{Symbol}, params::Vector{Symbol}) where E
        new{E}(expr, vars, params, Dict{UInt64, Tuple{Num, RuntimeGeneratedFunction}}())
    end
end

function ThermoFactory(expr, vars=[:T, :P, :t, :x, :y, :z])
    newvars, params = extract_vars_params(expr, vars)
    ThermoFactory{typeof(expr)}(expr, newvars, params)
end

ThermoFactory(sym::Symbol) = ThermoFactory{Symbol}(sym, [sym], Symbol[])

"""
    (factory::ThermoFactory)(; kwargs...)

Create a ThermoFunction with caching for optimal performance.
"""
function (factory::ThermoFactory)(; kwargs...)
    param_vals = Dict{Symbol, Any}(p => get(kwargs, p, 0.0) for p in factory.params)
    refs = NamedTuple([v => kwargs[v] for v in factory.vars if haskey(kwargs, v)])
    cache_key = hash(tuple(sort(collect(pairs(param_vals)))...))

    simplified, compiled = get!(factory.cache, cache_key) do
        var_sym_dict = Dict{Symbol, Num}(v => Symbolics.variable(v) for v in factory.vars)
        param_sym_dict = Dict{Symbol, Num}(p => Symbolics.variable(p) for p in factory.params)
        all_symbols = merge(var_sym_dict, param_sym_dict)

        symbolic_expr = Symbolics.parse_expr_to_symbolic(factory.expr, all_symbols)
        substitutions = Dict(param_sym_dict[p] => param_vals[p] for p in factory.params)
        substituted = Symbolics.substitute(symbolic_expr, substitutions)
        simplified = Symbolics.simplify(Symbolics.expand(substituted))
        compiled = compile_symbolic(simplified, factory.vars)

        (simplified, compiled)
    end

    return ThermoFunction(simplified, Tuple(factory.vars), refs, compiled)
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

# OPTIMIZED callable - specialized for 1, 2, 3+ variables
@inline function (tf::ThermoFunction{1})(; kwargs...)
    # Single variable - specialized, zero allocation
    v = tf.vars[1]
    val = haskey(kwargs, v) ? kwargs[v] : tf.refs[v]
    return tf.compiled(val)
end

@inline function (tf::ThermoFunction{2})(; kwargs...)
    # Two variables - specialized
    v1, v2 = tf.vars
    val1 = haskey(kwargs, v1) ? kwargs[v1] : get(tf.refs, v1, nothing)
    val2 = haskey(kwargs, v2) ? kwargs[v2] : get(tf.refs, v2, nothing)
    return tf.compiled(val1, val2)
end

@inline function (tf::ThermoFunction{3})(; kwargs...)
    # Three variables - specialized
    v1, v2, v3 = tf.vars
    val1 = haskey(kwargs, v1) ? kwargs[v1] : get(tf.refs, v1, nothing)
    val2 = haskey(kwargs, v2) ? kwargs[v2] : get(tf.refs, v2, nothing)
    val3 = haskey(kwargs, v3) ? kwargs[v3] : get(tf.refs, v3, nothing)
    return tf.compiled(val1, val2, val3)
end

@inline function (tf::ThermoFunction{N})(; kwargs...) where N
    # General case - fallback for 4+ variables
    if isempty(kwargs)
        # Fast path: no overrides, use refs directly
        var_values = ntuple(i -> tf.refs[tf.vars[i]], N)
    else
        # Merge path
        merged = merge(tf.refs, kwargs)
        var_values = ntuple(i -> merged[tf.vars[i]], N)
    end
    return tf.compiled(var_values...)
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

    return ThermoFunction(simplified, all_vars, refs, compiled)
end

"""
    combine_symbolic(op, tf::ThermoFunction, x::Number)

Combine ThermoFunction with scalar.
"""
function combine_symbolic(op, tf::ThermoFunction, x::Number)
    combined = op(tf.symbolic, x)
    simplified = Symbolics.simplify(Symbolics.expand(combined))
    compiled = compile_symbolic(simplified, collect(tf.vars))

    return ThermoFunction(simplified, tf.vars, tf.refs, compiled)
end

"""
    combine_symbolic(op, x::Number, tf::ThermoFunction)

Combine scalar with ThermoFunction.
"""
function combine_symbolic(op, x::Number, tf::ThermoFunction)
    combined = op(x, tf.symbolic)
    simplified = Symbolics.simplify(Symbolics.expand(combined))
    compiled = compile_symbolic(simplified, collect(tf.vars))

    return ThermoFunction(simplified, tf.vars, tf.refs, compiled)
end

"""
    apply_symbolic(op, tf::ThermoFunction)

Apply unary operation.
"""
function apply_symbolic(op, tf::ThermoFunction)
    combined = op(tf.symbolic)
    simplified = Symbolics.simplify(Symbolics.expand(combined))
    compiled = compile_symbolic(simplified, collect(tf.vars))

    return ThermoFunction(simplified, tf.vars, tf.refs, compiled)
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

# Display
function Base.show(io::IO, tf::ThermoFunction)
    print(io, tf.symbolic)
    if !isempty(tf.refs)
        print(io, " ◆ ", join(["$k=$v" for (k, v) in pairs(tf.refs)], ", "))
    end
end

function Base.show(io::IO, ::MIME"text/plain", tf::ThermoFunction)
    println(io, "ThermoFunction:")
    print(io, "  Expression: ")
    println(io, tf.symbolic)

    if !isempty(tf.refs)
        print(io, "  References: ")
        println(io, join(["$k=$v" for (k, v) in pairs(tf.refs)], ", "))
    end

    print(io, "  Variables: ")
    print(io, join(tf.vars, ", "))
end

function Base.show(io::IO, factory::ThermoFactory)
    print(io, factory.expr)
    if !isempty(factory.params)
        print(io, " ◆ params = ", join(factory.params, ", "))
    end
    if !isempty(factory.vars)
        print(io, " ◆ vars = ", join(factory.vars, ", "))
    end
end

function Base.show(io::IO, ::MIME"text/plain", factory::ThermoFactory)
    println(io, "ThermoFactory:")
    print(io, "  Expression: ")
    println(io, factory.expr)

    if !isempty(factory.params)
        print(io, "  Parameters: ")
        println(io, join(factory.params, ", "))
    end

    print(io, "  Variables: ")
    print(io, join(factory.vars, ", "))
end

# Helper functions
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

# Constructors
function ThermoFunction(sym::Symbol; kwargs...)
    factory = ThermoFactory(sym)
    return factory(; kwargs...)
end

function ThermoFunction(expr, vars=[:T, :P, :t, :x, :y, :z]; kwargs...)
    factory = ThermoFactory(expr, vars)
    return factory(; kwargs...)
end

function ThermoFunction(symexpr::Union{SymbolicUtils.BasicSymbolic, Num}; kwargs...)
    vars = Symbol.(get_variables(symexpr))
    refs = NamedTuple([v => kwargs[v] for v in vars if haskey(kwargs, v)])
    compiled = compile_symbolic(symexpr, collect(vars))
    return ThermoFunction(Num(symexpr), Tuple(vars), refs, compiled)
end
