"""
    Callable

Abstract base type for all callable thermodynamic functions.
"""
abstract type Callable end

"""
    ThermoFunction{U,F,R} <: Callable

Representation of a thermodynamic function with symbolic expression, variables, units, and reference conditions.

# Fields

  - `symexpr::Num`: Symbolic expression of the thermodynamic function
  - `vars::Dict{Symbol,Num}`: Dictionary of variables (symbols to symbolic expressions)
  - `unit::U`: Unit of the function's output
  - `func::F`: Compiled function for evaluation
  - `ref::R`: Reference conditions dictionary

# Examples

```jldoctest
julia> tf = ThermoFunction(:Cp, [:a₀ => 1.0, :a₁ => 2.0])
ThermoFunction with expression: a₀ + a₁*T
```
"""
struct ThermoFunction{U,F,R} <: Callable
    symexpr::Num
    vars::Dict{Symbol,Num}
    unit::U
    func::F
    ref::R
end

"""
    ThermoFunction(symexpr::Num, vars::Dict{Symbol,Num}, unit::U, ref::R) where {U,R}

Construct a ThermoFunction from symbolic expression, variables, unit, and reference.

# Arguments

  - `symexpr`: Symbolic expression
  - `vars`: Dictionary of variables
  - `unit`: Unit of the function
  - `ref`: Reference conditions

# Examples

```jldoctest
julia> tf = ThermoFunction(expr, vars, u"J/(mol*K)", Dict(:T => 298.15u"K"))
expr = Num(:a₀ + :a₁ * :T)
```
"""
function ThermoFunction(symexpr::Num, vars::Dict{Symbol,Num}, unit::U, ref::R) where {U,R}
    func = eval(build_function(symexpr, values(vars)...; expression=Val{false}))
    return ThermoFunction(symexpr, vars, unit, func, ref)
end

"""
    thermo_function_library

Dictionary of predefined thermodynamic function expressions.

# Contents

  - `:Cp`: Heat capacity expression
  - `:CpoverT`: Heat capacity divided by temperature
  - `:logKr`: Logarithm of equilibrium constant

# Examples

```jldoctest
julia> thermo_function_library[:Cp]
:(a₀ + a₁ * T + a₂ / T^2 + a₃ / √T + a₄ * T^2 + a₅ * T^3 + a₆ * T^4 + a₇ / T^3 + a₈ / T + a₉ * √T + a₁₀ * log(T))
```
"""
const thermo_function_library = Dict(
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
    :CpoverT => :(
        a₀ / T +
        a₁ +
        a₂ / T^3 +
        a₃ / T^(3 / 2) +
        a₄ * T +
        a₅ * T^2 +
        a₆ * T^3 +
        a₇ / T^4 +
        a₈ / T^2 +
        a₉ / √T +
        a₁₀ * log(T) / T
    ),
    :logKr => :(A₀ + A₁ * T + A₂ / T + A₃ * log(T) + A₄ / T^2 + A₅ * T^2 + A₆ * √T),
)

"""
    ThermoFunction(expr::Union{Symbol,Expr}, params=Pair[], vars=[:T, :P, :t]; ref=[:T => 298.15u"K", :P => 1u"bar", :t => 0u"s"])

Construct a ThermoFunction from an expression, parameters, variables, and reference conditions.

# Arguments

  - `expr`: Symbol or expression (or key from thermo_function_library)
  - `params`: Parameter values as Pairs
  - `vars`: Variable symbols (default: [:T, :P, :t])
  - `ref`: Reference conditions dictionary

# Examples

```jldoctest
julia> tf = ThermoFunction(:Cp, [:a₀ => 1.0, :a₁ => 2.0])
ThermoFunction with Cp expression

julia> tf = ThermoFunction(:(a₀ + a₁*T), [:a₀ => 1.0, :a₁ => 2.0])
Custom thermodynamic function
```
"""
function ThermoFunction(
    expr::Union{Symbol,Expr},
    params=Pair[],
    vars=[:T, :P, :t];
    ref=[:T => 298.15u"K", :P => 1u"bar", :t => 0u"s"],
)
    expr = get(thermo_function_library, expr, expr)
    symexpr = Num(parse_expr_to_symbolic(expr, @__MODULE__))
    varofexpr = Num.(get_variables(symexpr))
    dictallvars = Dict(Symbol(v) => v for v in varofexpr)
    dictparams = Dict(dictallvars[k] => v for (k, v) in params if haskey(dictallvars, k))
    dictref = Dict(dictallvars[k] => v for (k, v) in ref if haskey(dictallvars, k))
    vecvars = filter(x -> Symbol(x) ∈ vars, varofexpr)
    vecparams = filter(x -> Symbol(x) ∉ vars, varofexpr)
    veczeros = filter(x -> x ∉ keys(dictparams), vecparams)
    symexpr = substitute(symexpr, Dict(k => 0 for k in veczeros))
    nounitfunc = vec([
        f(var) => 1 for f in [
            log,
            log10,
            log2,
            sin,
            cos,
            tan,
            asin,
            acos,
            atan,
            sinh,
            cosh,
            tanh,
            asinh,
            acosh,
            atanh,
            exp,
        ],
        var in vecvars
    ])
    unit = if promote_type(typeof.(values(dictparams))...) <: Quantity
        Quantity(
            1,
            dimension(
                substitute(symexpr, [collect(dictparams); nounitfunc; collect(dictref)])
            ),
        )
    else
        1
    end
    symexpr = substitute(symexpr, [k => ustrip(v) for (k, v) in dictparams])
    dictvars = Dict{Symbol,Num}(zip(Symbol.(vecvars), vecvars))
    func = eval(build_function(symexpr, vecvars...; expression=Val{false}))
    return ThermoFunction(symexpr, dictvars, unit, func, Dict(ref))
end

"""
    (tf::ThermoFunction)(vars...)

Evaluate the thermodynamic function at given variable values.

# Arguments

  - `vars`: Variable values in the same order as tf.vars

# Returns

  - The function value (scalar)

# Examples

```jldoctest
julia> tf(300.0)  # Evaluate at T=300
tf = ThermoFunction(:Cp, [:a₀ => 1.0, :a₁ => 2.0])
```
"""
function (tf::ThermoFunction)(vars...)
    return tf.func(vars...)
end

"""
    (tf::ThermoFunction)(vars::Quantity...)

Evaluate the thermodynamic function at given variable values with units.

# Arguments

  - `vars`: Variable values with units

# Returns

  - The function value with appropriate units

# Examples

```jldoctest
julia> tf(300.0u"K")  # Returnss value with units
tf = ThermoFunction(:Cp, [:a₀ => 1.0u"J/(mol*K)", :a₁ => 2.0u"J/(mol*K^2)"], [:T => 298.15u"K"])
```
"""
function (tf::ThermoFunction)(vars::Quantity...)
    return tf.func(ustrip.(vars)...) * tf.unit
end

"""
    (tf::ThermoFunction)()

Evaluate the thermodynamic function at reference conditions.

# Returns

  - The function value at reference conditions

# Examples

```jldoctest
julia> tf()  # Evaluate at reference T=298.15
tf = ThermoFunction(:Cp, [:a₀ => 1.0, :a₁ => 2.0], ref=Dict(:T => 298.15))
```
"""
function (tf::ThermoFunction)()
    vars = [tf.ref[var] for var in keys(tf.vars)]
    if length(vars) > 0
        return tf(vars...)
    else
        return tf.func(0)
    end
end

"""
    +(tf::ThermoFunction, x::Number)

Add a number to a thermodynamic function.

# Arguments

  - `tf`: ThermoFunction
  - `x`: Number to add (must have compatible units)

# Returns

  - New ThermoFunction with added constant

# Examples

```jldoctest
julia> tf2 = tf + 5.0
tf = ThermoFunction(:Cp, [:a₀ => 1.0])
```
"""
function +(tf::ThermoFunction, x::Number)
    @assert dimension(tf.unit) == dimension(x)
    return ThermoFunction(tf.symexpr + ustrip(x), tf.vars, tf.unit, tf.ref)
end

"""
    +(x::Number, tf::ThermoFunction)

Add a number to a thermodynamic function (commutative version).

# Examples

```jldoctest
julia> tf2 = 5.0 + tf
tf = ThermoFunction(:Cp, [:a₀ => 1.0])
```
"""
+(x::Number, tf::ThermoFunction) = +(tf, x)

"""
    -(tf::ThermoFunction)

Negate a thermodynamic function.

# Returns

  - New ThermoFunction with negated expression

# Examples

```jldoctest
julia> tf2 = -tf
tf = ThermoFunction(:Cp, [:a₀ => 1.0])
```
"""
-(tf::ThermoFunction) = ThermoFunction(-1 * tf.symexpr, tf.vars, tf.unit, tf.ref)

"""
    *(tf::ThermoFunction, x::Number)

Multiply a thermodynamic function by a number.

# Arguments

  - `tf`: ThermoFunction
  - `x`: Scaling factor

# Returns

  - New ThermoFunction with scaled expression and units

# Examples

```jldoctest
julia> tf2 = tf * 2.0
tf = ThermoFunction(:Cp, [:a₀ => 1.0])
```
"""
function *(tf::ThermoFunction, x::Number)
    return ThermoFunction(
        ModelingToolkit.expand(tf.symexpr * ustrip(x)),
        tf.vars,
        tf.unit * dimension(x),
        tf.ref,
    )
end

"""
    *(x::Number, tf::ThermoFunction)

Multiply a thermodynamic function by a number (commutative version).
"""
*(x::Number, tf::ThermoFunction) = *(tf, x)

"""
    -(tf::ThermoFunction, x)

Subtract a number from a thermodynamic function.

# Examples

```jldoctest
julia> tf2 = tf - 5.0
tf = ThermoFunction(:Cp, [:a₀ => 1.0])
```
"""
-(tf::ThermoFunction, x) = +(tf, -1 * x)

"""
    -(x::Number, tf::ThermoFunction)

Subtract a thermodynamic function from a number.
"""
-(x::Number, tf::ThermoFunction) = +(-tf, x)

"""
    /(tf::ThermoFunction, x::Number)

Divide a thermodynamic function by a number.

# Returns

  - New ThermoFunction with divided expression and units

# Examples

```jldoctest
julia> tf2 = tf / 2.0
tf = ThermoFunction(:Cp, [:a₀ => 1.0])
```
"""
function /(tf::ThermoFunction, x::Number)
    return ThermoFunction(
        ModelingToolkit.expand(tf.symexpr / ustrip(x)),
        tf.vars,
        tf.unit / dimension(x),
        tf.ref,
    )
end

"""
    /(x::Number, tf::ThermoFunction)

Divide a number by a thermodynamic function.
"""
function /(x::Number, tf::ThermoFunction)
    return ThermoFunction(
        ModelingToolkit.expand(ustrip(x) / tf.symexpr),
        tf.vars,
        dimension(x) / tf.unit,
        tf.ref,
    )
end

"""
    ^(tf::ThermoFunction, x::Number)

Raise a thermodynamic function to a power.

# Arguments

  - `x`: Power (must be dimensionless)

# Examples

```jldoctest
julia> tf2 = tf^2
tf = ThermoFunction(:Cp, [:a₀ => 1.0])
```
"""
function ^(tf::ThermoFunction, x::Number)
    @assert dimension(1) == dimension(x)
    return ThermoFunction(tf.symexpr^x, tf.vars, tf.unit^x, tf.ref)
end

"""
    contained_dict(d1, d2)

Check if one dictionary is contained within another.

# Arguments

  - `d1`: First dictionary
  - `d2`: Second dictionary

# Returns

  - The containing dictionary if one contains the other, nothing otherwise

# Examples

```jldoctest
julia> contained_dict(d1, d2) == d2
d1 = Dict(:a => 1, :b => 2)
```
"""
function contained_dict(d1, d2)
    c1 = all(k -> haskey(d2, k) && isequal(d2[k], d1[k]), keys(d1))
    c2 = all(k -> haskey(d1, k) && isequal(d1[k], d2[k]), keys(d2))
    if c1
        return d2
    elseif c2
        return d1
    else
        return nothing
    end
end

"""
    +(tf1::ThermoFunction, tf2::ThermoFunction)

Add two thermodynamic functions.

# Arguments

  - `tf1`: First ThermoFunction
  - `tf2`: Second ThermoFunction

# Returns

  - New ThermoFunction with combined expressions

# Examples

```jldoctest
julia> tf3 = tf1 + tf2
tf1 = ThermoFunction(:Cp, [:a₀ => 1.0])
```
"""
function +(tf1::ThermoFunction, tf2::ThermoFunction)
    vars = contained_dict(tf1.vars, tf2.vars)
    @assert !isnothing(vars) && dimension(tf1.unit) == dimension(tf2.unit)
    return ThermoFunction(
        tf1.symexpr + tf2.symexpr, vars, tf1.unit, merge(tf1.ref, tf2.ref)
    )
end

"""
    *(tf1::ThermoFunction, tf2::ThermoFunction)

Multiply two thermodynamic functions.

# Examples

```jldoctest
julia> tf3 = tf1 * tf2
tf1 = ThermoFunction(:Cp, [:a₀ => 1.0])
```
"""
function *(tf1::ThermoFunction, tf2::ThermoFunction)
    vars = contained_dict(tf1.vars, tf2.vars)
    @assert !isnothing(vars)
    return ThermoFunction(
        tf1.symexpr * tf2.symexpr, vars, tf1.unit * tf2.unit, merge(tf1.ref, tf2.ref)
    )
end

"""
    /(tf1::ThermoFunction, tf2::ThermoFunction)

Divide two thermodynamic functions.

# Examples

```jldoctest
julia> tf3 = tf1 / tf2
tf1 = ThermoFunction(:Cp, [:a₀ => 1.0])
```
"""
function /(tf1::ThermoFunction, tf2::ThermoFunction)
    vars = contained_dict(tf1.vars, tf2.vars)
    @assert !isnothing(vars)
    return ThermoFunction(
        tf1.symexpr / tf2.symexpr, vars, tf1.unit / tf2.unit, merge(tf1.ref, tf2.ref)
    )
end

"""
    Base.inv(tf::ThermoFunction)

Compute the inverse of a thermodynamic function.

# Examples

```jldoctest
julia> tf2 = inv(tf)
tf = ThermoFunction(:Cp, [:a₀ => 1.0])
```
"""
Base.inv(tf::ThermoFunction) = 1 / tf

"""
    ∂(tf::ThermoFunction, var=collect(keys(tf.vars))[1])

Compute the partial derivative of a thermodynamic function.

# Arguments

  - `var`: Variable to differentiate with respect to (default: first variable)

# Returns

  - New ThermoFunction representing the derivative

# Examples

```jldoctest
julia> dtf = ∂(tf)  # Derivative with respect to T
tf = ThermoFunction(:Cp, [:a₀ => 1.0, :a₁ => 2.0])
```
"""
function ∂(tf::ThermoFunction, var=collect(keys(tf.vars))[1])
    return ThermoFunction(
        Num(
            ModelingToolkit.expand(
                expand_derivatives(ModelingToolkit.Differential(tf.vars[var])(tf.symexpr))
            ),
        ),
        tf.vars,
        tf.unit / dimension(tf.ref[var]),
        tf.ref,
    )
end

"""
    ∫(tf::ThermoFunction, var=collect(keys(tf.vars))[1])

Compute the integral of a thermodynamic function.

# Arguments

  - `var`: Variable to integrate with respect to (default: first variable)

# Returns

  - New ThermoFunction representing the integral

# Examples

```jldoctest
julia> itf = ∫(tf)  # Integral with respect to T
tf = ThermoFunction(:Cp, [:a₀ => 1.0, :a₁ => 2.0])
```
"""
function ∫(tf::ThermoFunction, var=collect(keys(tf.vars))[1])
    return ThermoFunction(
        Num(
            ModelingToolkit.expand(
                SymbolicNumericIntegration.integrate(
                    tf.symexpr, tf.vars[var]; symbolic=true, detailed=false
                ),
            ),
        ),
        tf.vars,
        tf.unit * dimension(tf.ref[var]),
        tf.ref,
    )
end

"""
    Base.show(io::IO, tf::ThermoFunction)

Display a thermodynamic function in a readable format.

# Arguments

  - `io`: Output stream
  - `tf`: ThermoFunction to display

# Examples

```jldoctest
julia> show(stdout, tf)
tf = ThermoFunction(:Cp, [:a₀ => 1.0, :a₁ => 2.0])
```
"""
function Base.show(io::IO, tf::ThermoFunction)
    print(io, tf.symexpr)
    print(io, " ♢ unit=[", dimension(tf.unit), "]")
    print(io, " ♢ ref=[")
    print(io, join(["$(k)=$(v)" for (k, v) in tf.ref], ", "), "]")
end

"""
    apply(func::Function, tf::ThermoFunction, args...; kwargs...)

Apply a function to a thermodynamic function's expression.

# Arguments

  - `func`: Function to apply
  - `tf`: ThermoFunction
  - `args...`: Additional arguments for func
  - `kwargs...`: Additional keyword arguments

# Returns

  - New ThermoFunction with transformed expression

# Examples

```jldoctest
julia> tf2 = apply(exp, tf)
tf = ThermoFunction(:Cp, [:a₀ => 1.0, :a₁ => 2.0])
```
"""
function apply(func::Function, tf::ThermoFunction, args...; kwargs...)
    return ThermoFunction(
        Num(func(tf.symexpr, args...; kwargs...)), tf.vars, func(tf.unit), tf.ref
    )
end

const MATH_FUNCTIONS = [
    :expand,
    :simplify,
    :sqrt,
    :log,
    :log10,
    :log2,
    :sin,
    :cos,
    :tan,
    :asin,
    :acos,
    :atan,
    :sinh,
    :cosh,
    :tanh,
    :asinh,
    :acosh,
    :atanh,
    :exp,
]

for fn in MATH_FUNCTIONS
    mod = (fn in [:expand, :simplify]) ? ModelingToolkit : Base
    func_name = fn

    @eval begin
        function $(mod).$(func_name)(tf::ThermoFunction)
            return ThermoFunction(
                Num($(mod).$(func_name)(tf.symexpr)),
                tf.vars,
                $(mod).$(func_name)(tf.unit),
                tf.ref,
            )
        end
    end

    @eval @doc """
        $($fn)(tf::ThermoFunction)
    
    Apply the $($fn) function to a thermodynamic function's expression.
    
    # Return
    - New ThermoFunction with $($fn) applied to the expression and unit
    
    # Examples
    ```jldoctest
    julia> tf = ThermoFunction(:Cp, [:a0 => 1.0, :a1 => 2.0])
    julia> tf2 = $($fn)(tf)
    # Returns a new ThermoFunction with $($fn) applied to the expression
    ```
    """ $(mod).$(func_name)
end
