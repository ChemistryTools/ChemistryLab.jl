*(x::Num, y::AbstractQuantity) = Quantity(x * y.value, y.dimensions)
*(x::AbstractQuantity, y::Num) = *(y, x)

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
:(a₀ + a₁ * T + a₂ / T ^ 2 + a₃ / √T + a₄ * T ^ 2 + a₅ * T ^ 3 + a₆ * T ^ 4 + a₇ / T ^ 3 + a₈ / T + a₉ * √T + a₁₀ * log(T))

julia> thermo_function_library[:CpoverT]
:(a₀ / T + a₁ + a₂ / T ^ 3 + a₃ / T ^ (3 / 2) + a₄ * T + a₅ * T ^ 2 + a₆ * T ^ 3 + a₇ / T ^ 4 + a₈ / T ^ 2 + a₉ / √T + (a₁₀ * log(T)) / T)

julia> thermo_function_library[:logKr]
:(A₀ + A₁ * T + A₂ / T + A₃ * log(T) + A₄ / T ^ 2 + A₅ * T ^ 2 + A₆ * √T)
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
    ThermoFunction(expr::Union{Symbol,Expr}, params=Pair[], vars=[:T, :P, :t]; ref=[])

Construct a ThermoFunction from an expression, parameters, variables, and reference conditions.

# Arguments

  - `expr`: Symbol or expression (or key from thermo_function_library)
  - `params`: Parameter values as Pairs
  - `vars`: Variable symbols (default: [:T, :P, :t])
  - `ref`: Reference conditions dictionary (default: ref=[:T => 298.15u"K", :P => 1u"bar", :t => 0u"s"])

!!! warning
    It is possible to use or not units when defining parameters and reference values but in case of dimensioned quantities, the user should take care of the homogeneity of the expression consistently with parameters and variables units. However variables inside functions such as `log`, `exp`, `cos`... are assumed to be dimensionless (or better implicitly divided by a unit value in the International System of Units) so that it is possible to write `log(T)`, `exp(T)`, `cos(T)` while keeping `T` in Kelvins.

# Examples

```jldoctest
julia> expr = :(α + β * T + γ * log(T))
:(α + β * T + γ * log(T))

julia> params = [:α => 210.0u"J/K/mol", :β => 0.0u"J/mol/K^2", :γ => -3.07e6u"J/K/mol"]
3-element Vector{Pair{Symbol, Quantity{Float64, Dimensions{FRInt32}}}}:
 :α => 210.0 m² kg s⁻² K⁻¹ mol⁻¹
 :β => 0.0 m² kg s⁻² K⁻² mol⁻¹
 :γ => -3.07e6 m² kg s⁻² K⁻¹ mol⁻¹

julia> vars = [:T]
1-element Vector{Symbol}:
 :T

julia> tf = ThermoFunction(expr, params, vars; ref=[:T => 298.15u"K"])
210.0 - 3.07e6log(T) ♢ unit=[m² kg s⁻² K⁻¹ mol⁻¹] ♢ ref=[T=298.15 K]

julia> tf() # default evaluation at reference variable here Tref=298.15K
-1.7491411916797183e7 m² kg s⁻² K⁻¹ mol⁻¹

julia> tf(300.0u"K")
-1.7510402197194535e7 m² kg s⁻² K⁻¹ mol⁻¹

julia> tf(300.0)
-1.7510402197194535e7

julia> tf_nounits = ThermoFunction(expr, [:α => 210.0, :β => 0.0, :γ => -3.07e6], vars; ref=[:T => 298.15])
210.0 - 3.07e6log(T) ♢ unit=[] ♢ ref=[T=298.15]

julia> tf_nounits(300.0u"K")
-1.7510402197194535e7

julia> tf_nounits(300.)
-1.7510402197194535e7
```

!!! note
    Note that `ThermoFunction` can be used as a functor i.e. can apply on a dimensioned or dimensionless value, which respectively returns a dimensioned or dimensionless result. The default argument(s) is (are) the reference one(s).
"""
function ThermoFunction(
    expr::Union{Symbol,Expr},
    params=Pair[],
    vars=[:T, :P, :t];
    ref=[],
)
    expr = get(thermo_function_library, expr, expr)
    symexpr = Num(parse_expr_to_symbolic(expr, @__MODULE__))
    varofexpr = Num.(get_variables(symexpr))
    dictallvars = Dict(Symbol(v) => v for v in varofexpr)
    dictparams = Dict(dictallvars[k] => v for (k, v) in params if haskey(dictallvars, k))
    with_units = promote_type(typeof.(values(dictparams))...) <: Quantity
    default_ref = with_units ? Dict(:T => 298.15u"K", :P => 1u"bar", :t => 0u"s") :
                    Dict(:T => 298.15, :P => 100000., :t => 0.)
    givenref = Dict(ref)
    vecvars = filter(x -> Symbol(x) ∈ vars, varofexpr)
    dictvars = Dict{Symbol,Num}(zip(Symbol.(vecvars), vecvars))
    dictref = Dict(v=>get(givenref, v, get(default_ref, v, 0)) for v in union(keys(dictvars), keys(givenref)))
    dictvarref = Dict(dictvars[v]=>get(givenref, v, get(default_ref, v, 0)) for v in keys(dictvars))
    vecparams = filter(x -> Symbol(x) ∉ vars, varofexpr)
    veczeros = filter(x -> x ∉ keys(dictparams) || iszero(dictparams[x]), vecparams)
    if length(veczeros)==length(vecparams)
        veczeros = filter(x -> !isequal(x,dictallvars[first(params[1])]), vecparams)
    end
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
    unit = if with_units
        Quantity(
        1,
        dimension(substitute(symexpr, [nounitfunc; collect(dictparams); collect(dictvarref)])),
    )
    else
        1
    end
    symexpr = substitute(symexpr, [k => ustrip(v) for (k, v) in dictparams])
    func = eval(build_function(symexpr, vecvars...; expression=Val{false}))
    return ThermoFunction(symexpr, dictvars, unit, func, dictref)
end

"""
    (tf::ThermoFunction)(vars...)

Evaluate the thermodynamic function at given variable values.

# Arguments

  - `vars`: Variable values in the same order as tf.vars

# Returns

  - The function value (scalar)
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
"""
function (tf::ThermoFunction)(vars::Quantity...)
    return tf.func(ustrip.(vars)...) * tf.unit
end

"""
    (tf::ThermoFunction)()

Evaluate the thermodynamic function at reference conditions.

# Returns

  - The function value at reference conditions
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
"""
function +(tf::ThermoFunction, x::Number)
    @assert dimension(tf.unit) == dimension(x)
    return ThermoFunction(tf.symexpr + ustrip(x), tf.vars, tf.unit, tf.ref)
end

"""
    +(x::Number, tf::ThermoFunction)

Add a number to a thermodynamic function (commutative version).
"""
+(x::Number, tf::ThermoFunction) = +(tf, x)

"""
    -(tf::ThermoFunction)

Negate a thermodynamic function.

# Returns

  - New ThermoFunction with negated expression
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
julia> Cp = ThermoFunction(:Cp, [:a₀ => 210.0u"J/K/mol", :a₁ => 0.12u"J/mol/K^2", :a₂ => -3.07e6u"J*K/mol", :a₃ => 0.0u"J/mol/√K"]; ref=[:T => 298.15u"K", :P => 1u"bar"])
210.0 + 0.12T + -3.07e6 / (T^2) ♢ unit=[m² kg s⁻² K⁻¹ mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]

julia> ∂(Cp)
0.12 + 6.14e6 / (T^3) ♢ unit=[m² kg s⁻² K⁻² mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]

julia> rate = ThermoFunction(:((c₁+c₂*t)/(c₃+c₄*√t)), [:c₁ => 1.0, :c₂ => 2.0u"1/s", :c₃ => 3.0, :c₄ => 4.0u"1/√s"])
(1.0 + 2.0t) / (3.0 + 4.0sqrt(t)) ♢ unit=[] ♢ ref=[t=0.0 s]

julia> ∂(rate)(1u"h")
0.004165467097946903 s⁻¹
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
        oneunit(tf.unit / Quantity(1, dimension(tf.ref[var]))),
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
julia> Cp = ThermoFunction(:Cp, [:a₀ => 210.0u"J/K/mol", :a₁ => 0.12u"J/mol/K^2", :a₂ => -3.07e6u"J*K/mol", :a₃ => 0.0u"J/mol/√K"]; ref=[:T => 298.15u"K", :P => 1u"bar"])
210.0 + 0.12T + -3.07e6 / (T^2) ♢ unit=[m² kg s⁻² K⁻¹ mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]

julia> ∫Cp = ∫(Cp)
210.0T + 3.07e6 / T + 0.06(T^2) ♢ unit=[m² kg s⁻² mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]

julia> ΔfH⁰_Tref = -2.72e6u"J/mol"
-2.72e6 m² kg s⁻² mol⁻¹

julia> ΔfH⁰ = (ΔfH⁰_Tref -∫Cp()) + ∫Cp
-2.798241935804469e6 + 210.0T + 3.07e6 / T + 0.06(T^2) ♢ unit=[m² kg s⁻² mol⁻¹] ♢ ref=[T=298.15 K, P=100000.0 m⁻¹ kg s⁻²]
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
        oneunit(tf.unit * Quantity(1, dimension(tf.ref[var]))),
        tf.ref,
    )
end

"""
    Base.show(io::IO, tf::ThermoFunction)

Display a thermodynamic function in a readable format.

# Arguments

  - `io`: Output stream
  - `tf`: ThermoFunction to display
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

```julia
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

    # Returns
    - New ThermoFunction with $($fn) applied to the expression and unit

    # Examples
    ```julia
    julia> tf = ThermoFunction(:Cp, [:a0 => 1.0, :a1 => 2.0])
    julia> tf2 = $($fn)(tf)
    # Returns a new ThermoFunction with $($fn) applied to the expression
    ```
    """
    $(mod).$(func_name)
end

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
