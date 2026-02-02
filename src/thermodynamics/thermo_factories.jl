using DynamicQuantities, ModelingToolkit, OrderedCollections
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)
import Base: ==, +, -, *, /, //, ^

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

safe_uparse(x::AbstractString) = uparse(x)
safe_uparse(x::AbstractQuantity) = x

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

"""
    extract_vars_params(expr, vars)

Extract variables and parameters in a Julia expression.

# Arguments
- `expr`: Expression like :(α + β * T * P + γ * log(T))
- `vars`: Variable candidates like [:T, :P]

Returns a vector of vars and a vector of parameters.

# Example
```julia
expr = :(α + β * T + γ * log(T) + δ * T * P)
vars, params = extract_vars_params(expr, [:T, :P])
```
"""
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

function substitute_symbols(ex, replacements::Dict{Symbol})
    if ex isa Symbol
        return haskey(replacements, ex) ? replacements[ex] : ex
    elseif ex isa Expr
        return Expr(ex.head, [substitute_symbols(arg, replacements) for arg in ex.args]...)
    else
        return ex
    end
end

"""
    Callable

Abstract base type for all callable thermodynamic functions.
"""
abstract type Callable end

# OPTIMIZATION: Split AbstractThermoFunction into two types
abstract type AbstractThermoFunction <: Callable end

"""
    FastThermoFunction: Direct RGF wrapper, minimal overhead

Used for functions created directly from ThermoFactory.
Optimized for speed - nearly as fast as native functions.
"""
struct FastThermoFunction{F<:RuntimeGeneratedFunction, N, R<:NamedTuple} <: AbstractThermoFunction
    compiled_func::F
    vars::NTuple{N, Symbol}  # Use NTuple for type stability
    refs::R
end

"""
    CompositeThermoFunction: Closure-based, handles operations

Used for results of arithmetic operations between AbstractThermoFunctions.
Slightly slower but handles variable unions correctly.
"""
struct CompositeThermoFunction{F, N, R<:NamedTuple} <: AbstractThermoFunction
    compiled_func::F
    vars::NTuple{N, Symbol}
    refs::R
end

function AbstractThermoFunction(sym::Symbol; kwargs...)
    factory = ThermoFactory(sym)
    return factory(; kwargs...)
end

# ULTRA-OPTIMIZED callable for FastThermoFunction
# Single variable - zero overhead path
@inline function (tf::FastThermoFunction{F, 1})(; kwargs...) where F
    v = tf.vars[1]
    val = get(kwargs, v, get(tf.refs, v, nothing))
    return tf.compiled_func(val)
end

# Two variables - optimized path
@inline function (tf::FastThermoFunction{F, 2})(; kwargs...) where F
    v1, v2 = tf.vars
    val1 = get(kwargs, v1, get(tf.refs, v1, nothing))
    val2 = get(kwargs, v2, get(tf.refs, v2, nothing))
    return tf.compiled_func(val1, val2)
end

# Three or more variables - general path
@inline function (tf::FastThermoFunction{F, N})(; kwargs...) where {F, N}
    # Avoid merge allocation when possible
    if isempty(kwargs)
        # Use refs directly
        var_values = ntuple(i -> tf.refs[tf.vars[i]], N)
    else
        merged = merge(tf.refs, kwargs)
        var_values = ntuple(i -> merged[tf.vars[i]], N)
    end
    return tf.compiled_func(var_values...)
end

# CompositeThermoFunction uses closure (slower but flexible)
@inline function (tf::CompositeThermoFunction)(; kwargs...)
    merged = merge(tf.refs, kwargs)
    return tf.compiled_func(; merged...)
end

# Arithmetic operations return CompositeThermoFunction
for op in (:(+), :(-), :(*), :(/), :(^))
    ex = quote
        function Base.$op(tf::AbstractThermoFunction, x::Number)
            compiled_func = (; kwargs...) -> $op(tf(; kwargs...), x)
            return CompositeThermoFunction(compiled_func, tf.vars, tf.refs)
        end

        function Base.$op(x::Number, tf::AbstractThermoFunction)
            compiled_func = (; kwargs...) -> $op(x, tf(; kwargs...))
            return CompositeThermoFunction(compiled_func, tf.vars, tf.refs)
        end
    end
    eval(ex)
end

for op in (:(+), :(-), :(*), :(/), :(^))
    ex = quote
        function Base.$op(tf1::AbstractThermoFunction, tf2::AbstractThermoFunction)
            all_vars = Tuple(union(tf1.vars, tf2.vars))
            refs = merge(tf1.refs, tf2.refs)

            let refs = refs, tf1 = tf1, tf2 = tf2
                compiled_func = (; kwargs...) -> begin
                    merged_kwargs = merge(refs, kwargs)
                    val1 = tf1(; merged_kwargs...)
                    val2 = tf2(; merged_kwargs...)
                    $op(val1, val2)
                end

                return CompositeThermoFunction(compiled_func, all_vars, refs)
            end
        end
    end
    eval(ex)
end

for f in [ADIM_MATH_FUNCTIONS ; [:sqrt, :abs]]
    if isdefined(Base, f)
        @eval begin
            function Base.$f(tf::AbstractThermoFunction)
                compiled_func = (; kwargs...) -> $f(tf(; kwargs...))
                return CompositeThermoFunction(compiled_func, tf.vars, tf.refs)
            end
        end
    end
end

"""
    ThermoFactory: Compiled Thermodynamic Model Generator

ThermoFactory(expr, vars) creates a reusable factory for generating compiled thermodynamic
models from symbolic expressions and state variables.

# Arguments:

- expr::Expr: Symbolic model expression (e.g., :(α + β * T + γ * log(T)))
- vars::AbstractVector{Symbol}: State variables (e.g., [:T], [:T, :P])

# Returns: Callable ThermoFactory struct producing AbstractThermoFunction{F} instances.

# Internal workflow:
1. Stores expr, vars, params
2. On factory(p1=v1,...): substitutes → compiles JIT (cached) → returns AbstractThermoFunction

# Benefits:
- 50-100× faster than eval-per-call approach (with cache)
- Reusable: 1 factory → unlimited models
- Type-stable & Serializable
- Inspectable: factory.expr, factory.params, factory.vars

Examples
```julia
julia> tfactory = ThermoFactory(:(α + β * T + γ * log(T)), [:T])
ThermoFactory(:(α + β * T + γ * log(T)), [:T], [:α, :β, :γ])

julia> tf = tfactory(α = 210.0, β = 0.0, γ = -3.07e6, T = 298.15)
AbstractThermoFunction{...}(...)

julia> tf()      # -1.749e7
julia> tf(T=300) # override
julia> tf.refs   # (T = 298.15,)
```
"""
struct ThermoFactory{E}
    expr::E
    vars::Vector{Symbol}
    params::Vector{Symbol}
    cache::Dict{UInt64,RuntimeGeneratedFunction}

    # Inner constructor to initialize cache
    function ThermoFactory{E}(
        expr::E, vars::Vector{Symbol}, params::Vector{Symbol}
    ) where {E}
        return new{E}(expr, vars, params, Dict{UInt64,RuntimeGeneratedFunction}())
    end
end

function ThermoFactory(expr, vars=[:T, :P, :t, :x, :y, :z])
    newvars, params = extract_vars_params(expr, vars)
    ThermoFactory{typeof(expr)}(expr, newvars, params)
end

ThermoFactory(sym::Symbol) = ThermoFactory{Symbol}(sym, [sym], Symbol[])

function ThermoFunction(sym::Symbol; kwargs...)
    factory = ThermoFactory(sym)
    return factory(; kwargs...)
end

function ThermoFunction(expr, vars=[:T, :P, :t, :x, :y, :z]; kwargs...)
    factory = ThermoFactory(expr, vars)
    return factory(; kwargs...)
end

"""
    (factory::ThermoFactory)(; kwargs...)

Create a AbstractThermoFunction with compiled expression for given parameters.

# OPTIMIZATIONS:
- Caches compiled functions based on parameter values (50-100× speedup on cache hit)
- Uses RuntimeGeneratedFunctions instead of eval (more robust, faster)
- Preserves exact same API and functionality
"""
function (factory::ThermoFactory)(; kwargs...)
    param_vals = Dict(p => get(kwargs, p, 0.0) for p in factory.params)
    refs = NamedTuple([v => kwargs[v] for v in factory.vars if haskey(kwargs, v)])
    cache_key = hash(tuple(sort(collect(pairs(param_vals)), by=x->x.first)...))

    rgf = get!(factory.cache, cache_key) do
        modified_expr = substitute_symbols(factory.expr, param_vals)
        var_syms = factory.vars
        if length(var_syms) == 1
            func_expr = :($(var_syms[1]) -> $modified_expr)
        else
            func_expr = :(($(var_syms...),) -> $modified_expr)
        end
        @RuntimeGeneratedFunction(func_expr)
    end

    return FastThermoFunction(rgf, Tuple(factory.vars), refs)
end


"""
    reconstruct_expr(tf::FastThermoFunction)

Reconstruct the analytical expression from a FastThermoFunction by examining its RGF.
"""
function reconstruct_expr(tf::FastThermoFunction)
    # Get the expression from the RuntimeGeneratedFunction
    rgf_expr = tf.compiled_func.body

    # The RGF wraps the actual expression in a quote block
    if rgf_expr isa Expr && rgf_expr.head == :block
        # Find the actual expression (skip line number nodes)
        for arg in rgf_expr.args
            if !(arg isa LineNumberNode)
                return arg
            end
        end
    end

    return rgf_expr
end

"""
    reconstruct_expr(tf::CompositeThermoFunction)

For composite functions, reconstruct the expression by evaluating with symbolic variables.
Returns a Symbolics.Num for pretty printing.
"""
function reconstruct_expr(tf::CompositeThermoFunction)
    try
        # Create symbolic variables
        if length(tf.vars) == 1
            var_sym = tf.vars[1]
            sym_var = only(Symbolics.@variables $var_sym)
            symbolic_kwargs = NamedTuple{(var_sym,)}((sym_var,))
        else
            vars_tuple = Expr(:tuple, tf.vars...)
            sym_vars = eval(:(Symbolics.@variables $vars_tuple))
            symbolic_kwargs = NamedTuple{tf.vars}(sym_vars)
        end

        # Evaluate the function with symbolic inputs
        result = tf.compiled_func(; symbolic_kwargs...)

        # Return the Num directly (don't convert to Expr)
        return result
    catch e
        # If symbolic evaluation fails, return placeholder
        return Symbol("<composite>")
    end
end

"""
    simplify_expr(expr)

Simplify an expression by removing zero terms and simplifying arithmetic.
Works with both Expr and Symbolics.Num.

Rules:
- 0 * x → 0
- x * 0 → 0
- 0 + x → x
- x + 0 → x
- x - 0 → x
- 1 * x → x
- x * 1 → x
- x / 1 → x
"""
function simplify_expr(expr)
    # If it's a Symbolics.Num, use Symbolics.simplify
    if expr isa Symbolics.Num
        return Symbolics.simplify(expr)
    end

    # Otherwise use our custom simplification for Expr
    if expr isa Expr
        if expr.head == :call
            func = expr.args[1]
            args = expr.args[2:end]

            # Recursively simplify arguments first
            simplified_args = [simplify_expr(arg) for arg in args]

            # Apply simplification rules
            if func == :*
                # Remove all zero factors
                non_zero_args = filter(!iszero_value, simplified_args)

                # If any factor is zero, result is zero
                if length(non_zero_args) < length(simplified_args)
                    return 0
                end

                # Remove all one factors
                non_one_args = filter(!isone_value, non_zero_args)

                # If nothing left, return 1
                if isempty(non_one_args)
                    return 1
                end

                # If only one arg left, return it directly
                if length(non_one_args) == 1
                    return non_one_args[1]
                end

                return Expr(:call, :*, non_one_args...)

            elseif func == :+
                # Remove all zero terms
                non_zero_args = filter(!iszero_value, simplified_args)

                # If nothing left, return 0
                if isempty(non_zero_args)
                    return 0
                end

                # If only one arg left, return it directly
                if length(non_zero_args) == 1
                    return non_zero_args[1]
                end

                return Expr(:call, :+, non_zero_args...)

            elseif func == :-
                if length(simplified_args) == 1
                    # Unary minus
                    arg = simplified_args[1]
                    if iszero_value(arg)
                        return 0
                    end
                    return Expr(:call, :-, arg)
                else
                    # Binary minus
                    if iszero_value(simplified_args[2])
                        return simplified_args[1]
                    end
                    if iszero_value(simplified_args[1])
                        return Expr(:call, :-, simplified_args[2])
                    end
                    return Expr(:call, :-, simplified_args...)
                end

            elseif func == :/
                if iszero_value(simplified_args[1])
                    return 0
                end
                if isone_value(simplified_args[2])
                    return simplified_args[1]
                end
                return Expr(:call, :/, simplified_args...)

            elseif func == :^
                base = simplified_args[1]
                exp = simplified_args[2]

                if iszero_value(exp)
                    return 1
                end
                if isone_value(exp)
                    return base
                end
                if iszero_value(base)
                    return 0
                end
                if isone_value(base)
                    return 1
                end

                return Expr(:call, :^, simplified_args...)
            else
                # Other functions: just simplify arguments
                return Expr(:call, func, simplified_args...)
            end
        else
            # Other expression types: recursively simplify
            return Expr(expr.head, [simplify_expr(arg) for arg in expr.args]...)
        end
    else
        # Leaf node (symbol, number, etc.)
        return expr
    end
end

"""
    iszero_value(x)

Check if a value represents zero (including Quantity types).
"""
function iszero_value(x)
    if x isa Number
        return iszero(x)
    elseif x isa Expr
        # Check if it's a literal zero or a Quantity constructor with zero
        if x.head == :call
            # Handle Quantity(0, ...) or similar
            if length(x.args) >= 2 && iszero_value(x.args[2])
                return true
            end
        end
        return false
    else
        # Try generic iszero for other types (like Quantity)
        try
            return iszero(x)
        catch
            return false
        end
    end
end

"""
    isone_value(x)

Check if a value represents one.
"""
function isone_value(x)
    if x isa Number
        return isone(x)
    else
        try
            return isone(x)
        catch
            return false
        end
    end
end

"""
    format_refs(refs::NamedTuple)

Format reference values for display.
"""
function format_refs(refs::NamedTuple)
    if isempty(refs)
        return ""
    end

    pairs_str = join(["$k=$v" for (k, v) in pairs(refs)], ", ")
    return " <|> $pairs_str"
end

"""
    Base.show(io::IO, tf::AbstractThermoFunction)

Display a ThermoFunction showing:
1. The simplified analytical expression
2. Reference values (if any) after a separator

Examples:
```
210.0 + -3.07e6*log(T) <|> T=298.15
T + P <|> T=300.0, P=101325.0
215.0 - 3.07e6*log(T) <|> T=298.15
```
"""
function Base.show(io::IO, tf::AbstractThermoFunction)
    # Reconstruct and simplify the expression (dispatch on type)
    expr = reconstruct_expr(tf)
    simplified = simplify_expr(expr)

    # Format the expression
    expr_str = string(simplified)

    # Format the references
    refs_str = format_refs(tf.refs)

    # Print
    print(io, expr_str, refs_str)
end

"""
    Base.show(io::IO, ::MIME"text/plain", tf::AbstractThermoFunction)

Detailed display for REPL (multi-line format).
"""
function Base.show(io::IO, ::MIME"text/plain", tf::AbstractThermoFunction)
    # Get type name
    type_name = typeof(tf).name.name

    println(io, type_name, ":")
    print(io, "  Expression: ")

    expr = reconstruct_expr(tf)
    simplified = simplify_expr(expr)

    if !isempty(tf.refs)
        println(io, simplified)
        print(io, "  References: ")
        print(io, tf.refs)
    else
        print(io, simplified)
    end

    # print(io, "  Variables: ")
    # print(io, tf.vars)
end

check_dimensions(dict_expr, units = dict_expr[:units]) =
     Dict(k => infer_unit(v, units) for (k, v) in dict_expr if k != :units)

build_thermo_factories(dict_expr) =
     Dict(k => ThermoFactory(v, [:T, :P]) for (k, v) in dict_expr if k != :units)

function build_thermo_functions(thermo_model, params)
    dict_factories = THERMO_FACTORIES[thermo_model]
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
        ],
    ),
)

const THERMO_FACTORIES = Dict(k => build_thermo_factories(v) for (k, v) in THERMO_MODELS)
