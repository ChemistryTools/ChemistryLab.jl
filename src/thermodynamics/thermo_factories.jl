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

"""
    ThermoFunction

A compiled thermodynamic function with optimized performance.
All operations (arithmetic, mathematical functions) return new ThermoFunctions
with freshly compiled RuntimeGeneratedFunctions.
"""
struct ThermoFunction{F<:RuntimeGeneratedFunction, N, R<:NamedTuple} <: Callable
    compiled_func::F
    vars::NTuple{N, Symbol}
    refs::R
end

"""
    ThermoFactory: Compiled Thermodynamic Model Generator

ThermoFactory(expr, vars) creates a reusable factory for generating compiled thermodynamic
models from symbolic expressions and state variables.

# Arguments:

- expr::Expr: Symbolic model expression (e.g., :(α + β * T + γ * log(T)))
- vars::AbstractVector{Symbol}: State variables (e.g., [:T], [:T, :P])

# Returns: Callable ThermoFactory struct producing ThermoFunction instances.

# Internal workflow:
1. Stores expr, vars, params
2. On factory(p1=v1,...): substitutes → compiles JIT (cached) → returns ThermoFunction

# Benefits:
- 50-100× faster than eval-per-call approach (with cache)
- Reusable: 1 factory → unlimited models
- Type-stable & generic: works with Float64, Quantity, Symbolics.Num
- Inspectable: factory.expr, factory.params, factory.vars

# Examples
```julia
julia> tfactory = ThermoFactory(:(α + β * T + γ * log(T)), [:T])
julia> tf = tfactory(α = 210.0, β = 0.0, γ = -3.07e6, T = 298.15)
julia> tf()      # Evaluate at reference T=298.15
julia> tf(T=300) # Override reference value
```
"""
struct ThermoFactory{E}
    expr::E
    vars::Vector{Symbol}
    params::Vector{Symbol}
    cache::Dict{UInt64, RuntimeGeneratedFunction}

    function ThermoFactory{E}(expr::E, vars::Vector{Symbol}, params::Vector{Symbol}) where E
        new{E}(expr, vars, params, Dict{UInt64, RuntimeGeneratedFunction}())
    end
end

function ThermoFactory(expr, vars=[:T, :P, :t, :x, :y, :z])
    newvars, params = extract_vars_params(expr, vars)
    ThermoFactory{typeof(expr)}(expr, newvars, params)
end

ThermoFactory(sym::Symbol) = ThermoFactory{Symbol}(sym, [sym], Symbol[])

function (factory::ThermoFactory)(; kwargs...)
    param_vals = Dict{Symbol, Any}(p => get(kwargs, p, 0.0) for p in factory.params)
    refs = NamedTuple([v => kwargs[v] for v in factory.vars if haskey(kwargs, v)])
    cache_key = hash(tuple(sort(collect(pairs(param_vals)))...))

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

    return ThermoFunction(rgf, Tuple(factory.vars), refs)
end

function ThermoFunction(sym::Symbol; kwargs...)
    factory = ThermoFactory(sym)
    return factory(; kwargs...)
end

function ThermoFunction(expr, vars=[:T, :P, :t, :x, :y, :z]; kwargs...)
    factory = ThermoFactory(expr, vars)
    return factory(; kwargs...)
end

# Callable interface
@inline function (tf::ThermoFunction)(; kwargs...)
    merged = merge(tf.refs, kwargs)
    var_values = ntuple(i -> merged[tf.vars[i]], length(tf.vars))
    return tf.compiled_func(var_values...)
end

"""
    reconstruct_expr(tf::ThermoFunction)

Reconstruct the analytical expression from a ThermoFunction.
Works by examining the RuntimeGeneratedFunction body.
"""
function reconstruct_expr(tf::ThermoFunction)
    rgf_expr = tf.compiled_func.body

    if rgf_expr isa Expr && rgf_expr.head == :block
        for arg in rgf_expr.args
            if !(arg isa LineNumberNode)
                return arg
            end
        end
    end

    return rgf_expr
end

"""
    simplify_expr(expr)

Aggressively simplify expressions with multiple passes.
"""
function simplify_expr(expr)
    prev = nothing
    current = expr
    max_iterations = 10
    iteration = 0

    # Iterate until no more changes or max iterations
    while prev != current && iteration < max_iterations
        prev = current
        current = eval_constants(current)
        current = simplify_algebraic(current)
        current = normalize_negations(current)
        current = flatten_operations(current)
        iteration += 1
    end

    return current
end

"""
    eval_constants(expr)

Evaluate all numeric sub-expressions immediately.
"""
function eval_constants(expr)
    if expr isa Number
        return expr
    elseif expr isa Symbol
        return expr
    elseif expr isa Expr && expr.head == :call
        func = expr.args[1]
        args = [eval_constants(arg) for arg in expr.args[2:end]]

        # If all arguments are numbers, evaluate immediately
        if all(arg -> arg isa Number, args)
            try
                result = eval(Expr(:call, func, args...))
                # Clean up floating point errors
                if result isa AbstractFloat
                    if abs(result) < 1e-14
                        return 0.0
                    elseif abs(result - round(result)) < 1e-10
                        return round(result)
                    end
                end
                return result
            catch
                return Expr(:call, func, args...)
            end
        else
            return Expr(:call, func, args...)
        end
    else
        return expr
    end
end

"""
    normalize_negations(expr)

Normalize all negations: turn subtractions and negative multiplications into additions.
- (a - b) → (a + (-1)*b)
- (-1)*x → -x
- -(-x) → x
- a - (-b) → a + b
"""
function normalize_negations(expr)
    if expr isa Expr && expr.head == :call
        func = expr.args[1]
        args = [normalize_negations(arg) for arg in expr.args[2:end]]

        if func == :-
            if length(args) == 1
                # Unary minus: -x
                arg = args[1]

                # -(-x) → x
                if arg isa Expr && arg.head == :call && arg.args[1] == :- && length(arg.args) == 2
                    return arg.args[2]
                end

                # -(a * b) where a is -1 → a * b with -1 removed
                if arg isa Expr && arg.head == :call && arg.args[1] == :*
                    mul_args = arg.args[2:end]
                    # Count negatives
                    neg_count = count(a -> (a isa Number && a == -1) || (a == -1.0), mul_args)
                    if neg_count > 0
                        # Remove -1s and add one more negation
                        filtered = filter(a -> !((a isa Number && a == -1) || (a == -1.0)), mul_args)
                        if isempty(filtered)
                            return -1
                        elseif length(filtered) == 1
                            return Expr(:call, :-, filtered[1])
                        else
                            return Expr(:call, :-, Expr(:call, :*, filtered...))
                        end
                    end
                end

                # -number → evaluate
                if arg isa Number
                    return -arg
                end

                return Expr(:call, :-, arg)
            else
                # Binary minus: a - b → a + (-b)
                left = args[1]
                right = args[2]

                # a - (-b) → a + b
                if right isa Expr && right.head == :call && right.args[1] == :- && length(right.args) == 2
                    return normalize_negations(Expr(:call, :+, left, right.args[2]))
                end

                # Convert to addition
                neg_right = normalize_negations(Expr(:call, :-, right))
                return normalize_negations(Expr(:call, :+, left, neg_right))
            end

        elseif func == :*
            # Remove -1 factors and count them
            neg_count = 0
            result_args = []

            for arg in args
                if (arg isa Number && arg == -1) || (arg isa Number && arg == -1.0)
                    neg_count += 1
                elseif arg isa Number && arg < 0
                    neg_count += 1
                    push!(result_args, -arg)
                else
                    push!(result_args, arg)
                end
            end

            # Apply negations
            if isempty(result_args)
                return (-1)^neg_count
            end

            result = length(result_args) == 1 ? result_args[1] : Expr(:call, :*, result_args...)

            if neg_count % 2 == 1
                return Expr(:call, :-, result)
            else
                return result
            end
        else
            return Expr(:call, func, args...)
        end
    else
        return expr
    end
end

"""
    flatten_operations(expr)

Flatten nested + and * operations and simplify.
(a + b) + c → a + b + c
(a * b) * c → a * b * c
"""
function flatten_operations(expr)
    if expr isa Expr && expr.head == :call
        func = expr.args[1]
        args = [flatten_operations(arg) for arg in expr.args[2:end]]

        if func == :+
            # Flatten nested additions
            flattened = []
            for arg in args
                if arg isa Expr && arg.head == :call && arg.args[1] == :+
                    append!(flattened, arg.args[2:end])
                else
                    push!(flattened, arg)
                end
            end

            # Remove zeros
            flattened = filter(!iszero_value, flattened)

            # Combine numeric terms
            numbers = filter(x -> x isa Number, flattened)
            non_numbers = filter(x -> !(x isa Number), flattened)

            result_args = non_numbers
            if !isempty(numbers)
                sum_num = sum(numbers)
                if !iszero_value(sum_num)
                    push!(result_args, sum_num)
                end
            end

            if isempty(result_args)
                return 0
            elseif length(result_args) == 1
                return result_args[1]
            else
                return Expr(:call, :+, result_args...)
            end

        elseif func == :*
            # Flatten nested multiplications
            flattened = []
            for arg in args
                if arg isa Expr && arg.head == :call && arg.args[1] == :*
                    append!(flattened, arg.args[2:end])
                else
                    push!(flattened, arg)
                end
            end

            # Check for zeros
            if any(iszero_value, flattened)
                return 0
            end

            # Remove ones
            flattened = filter(!isone_value, flattened)

            # Combine numeric terms
            numbers = filter(x -> x isa Number, flattened)
            non_numbers = filter(x -> !(x isa Number), flattened)

            result_args = non_numbers
            if !isempty(numbers)
                prod_num = prod(numbers)
                if iszero_value(prod_num)
                    return 0
                elseif !isone_value(prod_num)
                    push!(result_args, prod_num)
                end
            end

            if isempty(result_args)
                return 1
            elseif length(result_args) == 1
                return result_args[1]
            else
                return Expr(:call, :*, result_args...)
            end

        elseif func == :^
            base = args[1]
            exp = args[2]

            if iszero_value(exp)
                return 1
            elseif isone_value(exp)
                return base
            elseif iszero_value(base)
                return 0
            elseif isone_value(base)
                return 1
            elseif base isa Number && exp isa Number
                return base ^ exp
            else
                return Expr(:call, :^, base, exp)
            end

        elseif func == :/
            num = args[1]
            den = args[2]

            if iszero_value(num)
                return 0
            elseif isone_value(den)
                return num
            elseif num isa Number && den isa Number
                return num / den
            else
                return Expr(:call, :/, num, den)
            end
        else
            return Expr(:call, func, args...)
        end
    else
        return expr
    end
end

"""
    simplify_algebraic(expr)

Apply algebraic identities.
"""
function simplify_algebraic(expr)
    if expr isa Expr && expr.head == :call
        func = expr.args[1]
        args = [simplify_algebraic(arg) for arg in expr.args[2:end]]

        # T * log(T) - T can sometimes be factored, etc.
        # For now, just recursively simplify
        return Expr(:call, func, args...)
    else
        return expr
    end
end

function iszero_value(x)
    if x isa Number
        return abs(x) < 1e-14
    else
        try
            return iszero(x)
        catch
            return false
        end
    end
end

function isone_value(x)
    if x isa Number
        return abs(x - 1) < 1e-14
    else
        try
            return isone(x)
        catch
            return false
        end
    end
end

"""
    combine_rgf(op::Symbol, tf::ThermoFunction)

Apply a unary operator or function to a ThermoFunction, returning a new
ThermoFunction with a compiled RGF for maximum performance.
"""
function combine_rgf(op::Symbol, tf::ThermoFunction)
    expr = reconstruct_expr(tf)
    simplified = simplify_expr(expr)
    combined_expr = Expr(:call, op, simplified)

    if length(tf.vars) == 1
        func_expr = :($(tf.vars[1]) -> $combined_expr)
    else
        func_expr = :(($(tf.vars...),) -> $combined_expr)
    end

    rgf = @RuntimeGeneratedFunction(func_expr)
    return ThermoFunction(rgf, tf.vars, tf.refs)
end

"""
    combine_rgf(op::Symbol, tf1::ThermoFunction, tf2::ThermoFunction)

Combine two ThermoFunctions with a binary operator, returning a new ThermoFunction.
"""
function combine_rgf(op::Symbol, tf1::ThermoFunction, tf2::ThermoFunction)
    all_vars = Tuple(union(tf1.vars, tf2.vars))
    refs = merge(tf1.refs, tf2.refs)

    expr1 = reconstruct_expr(tf1)
    expr2 = reconstruct_expr(tf2)
    simplified1 = simplify_expr(expr1)
    simplified2 = simplify_expr(expr2)
    combined_expr = Expr(:call, op, simplified1, simplified2)

    if length(all_vars) == 1
        func_expr = :($(all_vars[1]) -> $combined_expr)
    else
        func_expr = :(($(all_vars...),) -> $combined_expr)
    end

    rgf = @RuntimeGeneratedFunction(func_expr)
    return ThermoFunction(rgf, all_vars, refs)
end

"""
    combine_rgf(op::Symbol, tf::ThermoFunction, x::Number)

Combine a ThermoFunction with a scalar.
"""
function combine_rgf(op::Symbol, tf::ThermoFunction, x::Number)
    expr = reconstruct_expr(tf)
    simplified = simplify_expr(expr)
    combined_expr = Expr(:call, op, simplified, x)

    if length(tf.vars) == 1
        func_expr = :($(tf.vars[1]) -> $combined_expr)
    else
        func_expr = :(($(tf.vars...),) -> $combined_expr)
    end

    rgf = @RuntimeGeneratedFunction(func_expr)
    return ThermoFunction(rgf, tf.vars, tf.refs)
end

"""
    combine_rgf(op::Symbol, x::Number, tf::ThermoFunction)

Combine a scalar with a ThermoFunction.
"""
function combine_rgf(op::Symbol, x::Number, tf::ThermoFunction)
    expr = reconstruct_expr(tf)
    simplified = simplify_expr(expr)
    combined_expr = Expr(:call, op, x, simplified)

    if length(tf.vars) == 1
        func_expr = :($(tf.vars[1]) -> $combined_expr)
    else
        func_expr = :(($(tf.vars...),) -> $combined_expr)
    end

    rgf = @RuntimeGeneratedFunction(func_expr)
    return ThermoFunction(rgf, tf.vars, tf.refs)
end

# Binary operations between ThermoFunctions
for op in (:+, :-, :*, :/, :^)
    @eval Base.$op(tf1::ThermoFunction, tf2::ThermoFunction) = combine_rgf($(QuoteNode(op)), tf1, tf2)
end

# Binary operations with scalars
for op in (:+, :-, :*, :/, :^)
    @eval begin
        Base.$op(tf::ThermoFunction, x::Number) = combine_rgf($(QuoteNode(op)), tf, x)
        Base.$op(x::Number, tf::ThermoFunction) = combine_rgf($(QuoteNode(op)), x, tf)
    end
end

# Unary minus
Base.:-(tf::ThermoFunction) = combine_rgf(:-, tf)

# Mathematical functions
for f in [ADIM_MATH_FUNCTIONS ; :sqrt]
    if isdefined(Base, f)
        @eval Base.$f(tf::ThermoFunction) = combine_rgf($(QuoteNode(f)), tf)
    end
end

# Display functions
function format_refs(refs::NamedTuple)
    isempty(refs) ? "" : " ◆ " * join(["$k=$v" for (k, v) in pairs(refs)], ", ")
end

function Base.show(io::IO, tf::ThermoFunction)
    expr = reconstruct_expr(tf)
    simplified = simplify_expr(expr)
    print(io, string(simplified), format_refs(tf.refs))
end

function Base.show(io::IO, ::MIME"text/plain", tf::ThermoFunction)
    println(io, "ThermoFunction:")
    print(io, "  Expression: ")
    expr = reconstruct_expr(tf)
    simplified = simplify_expr(expr)
    println(io, simplified)

    if !isempty(tf.refs)
        print(io, "  References: ")
        println(io, tf.refs)
    end

    print(io, "  Variables: ")
    print(io, tf.vars)
end
