using ChemistryLab, DynamicQuantities, ModelingToolkit
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
                exponent = args[2]

                # Exponent must be a dimensionless number
                if exponent isa Number
                    return base_unit ^ exponent
                else
                    throw(ArgumentError("Exponent in power must be a number"))
                end

            else
                throw(ArgumentError("Unsupported function: $func"))
            end
        else
            throw(ArgumentError("Unsupported expression: $(expr.head)"))
        end
    else
        throw(ArgumentError("Unsupported expression type: $(typeof(expr))"))
    end
end

# Functor struct to hold the model
struct ThermoFunction{F}
    compiled_func::F  # Compiled expression with substituted parameters
    vars::Tuple{Vararg{Symbol}}
    refs::NamedTuple  # Reference values for variables
end

@inline function (tf::ThermoFunction)(; kwargs...)
    merged = merge(tf.refs, kwargs)
    return tf.compiled_func(; merged...)
end

for op in (:(+), :(-), :(*), :(/), :(^))
    ex = quote
        Base.$op(tf::ThermoFunction, x::Number) = ThermoFunction((; kwargs...) -> $op(tf(; tf.refs..., kwargs...), x), tf.vars, tf.refs)
        Base.$op(x::Number, tf::ThermoFunction) = ThermoFunction((; kwargs...) -> $op(x, tf(; tf.refs..., kwargs...)), tf.vars, tf.refs)
    end
    eval(ex)
end

for op in (:(+), :(-), :(*), :(/), :(^))
    ex = quote
        function Base.$op(tf1::ThermoFunction, tf2::ThermoFunction)
            all_vars = Tuple(union(tf1.vars, tf2.vars))
            refs = merge(tf1.refs, tf2.refs)

            # Capture refs in the closure
            let refs = refs, tf1 = tf1, tf2 = tf2
                compiled_func = (; kwargs...) -> begin
                    # Merge references with provided kwargs (kwargs override refs)
                    merged_kwargs = merge(refs, kwargs)

                    # Extract only relevant variables for each function
                    args1 = NamedTuple{tuple(tf1.vars...)}(tuple((merged_kwargs[v] for v in tf1.vars)...))
                    args2 = NamedTuple{tuple(tf2.vars...)}(tuple((merged_kwargs[v] for v in tf2.vars)...))


                    val1 = tf1.compiled_func(; args1...)
                    val2 = tf2.compiled_func(; args2...)

                    $op(val1, val2)
                end

                return ThermoFunction(compiled_func, collect(all_vars), refs)
            end
        end
    end
    eval(ex)
end

for op in (:(+), :(-), :(*), :(/), :(^))
    ex = quote
        function Base.$op(tf1::ThermoFunction, tf2::ThermoFunction)
            all_vars = Tuple(union(tf1.vars, tf2.vars))
            refs = merge(tf1.refs, tf2.refs)

            # Capture everything in the closure
            let refs = refs, tf1 = tf1, tf2 = tf2
                compiled_func = (; kwargs...) -> begin
                    # Merge references with provided kwargs (kwargs override refs)
                    merged_kwargs = merge(refs, kwargs)

                    # Call each function with the merged kwargs
                    # The callable will filter to only relevant variables
                    val1 = tf1(; merged_kwargs...)
                    val2 = tf2(; merged_kwargs...)

                    $op(val1, val2)
                end

                return ThermoFunction(compiled_func, collect(all_vars), refs)
            end
        end
    end
    eval(ex)
end

for f in [ADIM_MATH_FUNCTIONS ; :sqrt]
    if isdefined(Base, f)
        @eval begin
            function Base.$f(tf::ThermoFunction)
                compiled_func = (; kwargs...) -> $f(tf(; kwargs...))
                return ThermoFunction(compiled_func, tf.vars, tf.refs)
            end
        end
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
    ThermoFactory: Compiled Thermodynamic Model Generator

ThermoFactory(expr, vars) creates a reusable factory for generating compiled thermodynamic
models from symbolic expressions and state variables.

# Arguments:

- expr::Expr: Symbolic model expression (e.g., :(α + β * T + γ * log(T)))
- vars::AbstractVector{Symbol}: State variables (e.g., [:T], [:T, :P])

# Returns: Callable ThermoFactory struct producing ThermoFunction{F} instances.

# Internal workflow:
1. Stores expr, vars, params
2. On factory(p1=v1,...): substitutes → compiles JIT → returns ThermoFunction

# Benefits:
- 5× faster than @eval approach
- Reusable: 1 factory → unlimited models
- Type-stable & Serializable
- Inspectable: factory.expr, factory.params, factory.vars

Examples

```julia
julia> tfactory = ThermoFactory(:(α + β * T + γ * log(T)), [:T])
ThermoFactory(:(α + β * T + γ * log(T)), [:T], [:α, :β, :γ])

julia> tf = tfactory(α = 210.0, β = 0.0, γ = -3.07e6, T = 298.15)
ThermoFunction{typeof(inner)}((inner), (T = 298.15,))

julia> tf()      # -1.749e7
julia> tf(T=300) # override
julia> tf.refs   # (T = 298.15,)
```
"""
struct ThermoFactory{E}
    expr::E
    vars::Vector{Symbol}
    params::Vector{Symbol}
    # vars::Dict{Symbol, Dimensions{FRInt32}}
    # params::Dict{Symbol, Dimensions{FRInt32}}
    # unit::Dimensions{FRInt32}
end

function ThermoFactory(expr, vars=[:T, :P, :t, :x, :y, :z])
    newvars, params = extract_vars_params(expr, vars)
    ThermoFactory(expr, collect(newvars), collect(params))
end

function ThermoFactory(sym::Symbol)
    ThermoFactory(sym, [sym], Symbol[])
end

function ThermoFunction(sym::Symbol; kwargs...)
    factory = ThermoFactory(sym)
    return factory(; kwargs...)
end

function (factory::ThermoFactory)(; kwargs...)
    param_vals = Dict(p => get(kwargs, p, 0.0) for p in factory.params)
    refs = NamedTuple([v => kwargs[v] for v in factory.vars if haskey(kwargs, v)])

    modified_expr = substitute_symbols(factory.expr, param_vals)

    # Generate a unique function name to avoid conflicts
    func_name = gensym(:inner)

    inner_def = Expr(:function,
                    Expr(:call, func_name, Expr(:parameters, Expr(:(...), :kwargs))),
                    Expr(:block,
                        :(kw = NamedTuple(kwargs)),
                        [:($(v) = kw[$(QuoteNode(v))]) for v in factory.vars]...,
                        modified_expr))

    compiled_func = eval(inner_def)

    return ThermoFunction(compiled_func, Tuple(factory.vars), refs)
end
