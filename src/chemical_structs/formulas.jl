"""
    struct AtomGroup{T}

Simple container pairing an atomic symbol with a numeric coefficient.

# Fields

  - `coef::T`: numeric coefficient.
  - `sym::Symbol`: atomic symbol.

# Examples

```jldoctest
julia> AtomGroup(:H, 2)
AtomGroup{Int64}(2, :H)

julia> AtomGroup(:Ca)
AtomGroup{Int64}(1, :Ca)
```
"""
struct AtomGroup{T<:Number}
    coef::T
    sym::Symbol
end

"""
    Base.convert(::Type{AtomGroup}, sym::Symbol) -> AtomGroup

Convert a `Symbol` to a unit `AtomGroup` (coefficient = 1).

# Examples

```jldoctest
julia> convert(AtomGroup, :O)
AtomGroup{Int64}(1, :O)
```
"""
Base.convert(::Type{AtomGroup}, sym::Symbol) = AtomGroup(1, sym)

"""
    AtomGroup(sym::Symbol) -> AtomGroup
    AtomGroup(sym::Symbol, coef::T) where {T<:Number} -> AtomGroup{T}

Constructors for `AtomGroup`.

# Arguments

  - `sym`: atomic symbol.
  - `coef`: numeric coefficient (default 1).

# Examples

```jldoctest
julia> AtomGroup(:C)
AtomGroup{Int64}(1, :C)

julia> AtomGroup(:C, 3)
AtomGroup{Int64}(3, :C)
```
"""
AtomGroup(sym::Symbol) = AtomGroup(1, sym)
AtomGroup(sym::Symbol, coef::T) where {T<:Number} = AtomGroup(coef, sym)

"""
    struct Formula{T}

Canonical container for a chemical formula.

# Fields

  - `expr::String`: original input expression.
  - `phreeqc::String`: PHREEQC-compatible representation.
  - `unicode::String`: Unicode pretty representation.
  - `colored::String`: colored terminal representation.
  - `composition::OrderedDict{Symbol,T}`: mapping element symbol to coefficient.
  - `charge::Int8`: formal integer charge.

# Examples

```jldoctest
julia> f = Formula("H2O");

julia> composition(f)[:H]
2
```
"""
struct Formula{T<:Number}
    expr::String
    phreeqc::String
    unicode::String
    colored::String
    composition::OrderedDict{Symbol,T}
    charge::Int8
end

"""
    stoichtype(f::Formula{T}) where {T} -> Type{T}

Return the numeric stoichiometric coefficient type `T` for formula `f`.

# Examples

```jldoctest
julia> stoichtype(Formula("H2O"))
Int64
```
"""
stoichtype(f::Formula{T}) where {T} = T

"""
    expr(f::Formula) -> String

Return the original expression string stored in `f`.

# Examples

```jldoctest
julia> expr(Formula("H2O"))
"H2O"
```
"""
expr(f::Formula) = f.expr

"""
    phreeqc(f::Formula) -> String

Return the PHREEQC-compatible representation of `f`.

# Examples

```jldoctest
julia> phreeqc(Formula("H2O"))
"H2O"
```
"""
phreeqc(f::Formula) = f.phreeqc

"""
    unicode(f::Formula) -> String

Return the Unicode pretty representation of `f`.

# Examples

```jldoctest
julia> phreeqc(Formula("C3A"))
"C3A"
```
"""
unicode(f::Formula) = f.unicode

"""
    colored(f::Formula) -> String

Return the colored terminal representation of `f`.

# Examples

```@example
julia> colored(Formula("Ca(HSiO3)+"))
```
"""
colored(f::Formula) = f.colored

"""
    composition(f::Formula) -> OrderedDict{Symbol,T}

Return the composition mapping (element symbol => coefficient).

# Examples

```@example
julia> composition(Formula("Ca(HSiO3)+"))
```
"""
composition(f::Formula) = f.composition

"""
    charge(f::Formula) -> Int8

Return the formal integer charge of the formula.

# Examples

```jldoctest
julia> charge(Formula("Ca(HSiO3)+"))
1
```
"""
charge(f::Formula) = f.charge

"""
    Formula(expr::AbstractString="") -> Formula{T}

Parse an input chemical formula string and return a `Formula`. Supports simple formulas,
parentheses, hydrates, and common charge notations. Special tokens like `Zz` (charge placeholder)
and `e` (electron) are handled.

# Arguments

  - `expr`: input formula string.

# Examples

```jldoctest
julia> f = Formula("SO4-2");

julia> composition(f)[:S]
1
```
"""
function Formula(expr::AbstractString="")
    if expr == "Zz" || expr == "Zz+" || expr == "Zz⁺"
        composition = OrderedDict{Symbol,Int}()
        charge = 1
    elseif expr == "e" || expr == "e-" || expr == "e⁻"
        composition = OrderedDict{Symbol,Int}()
        charge = -1
    else
        composition = parse_formula(expr)
        charge = extract_charge(expr)
    end
    phreeqc_expr = unicode_to_phreeqc(expr)
    unicode_expr = phreeqc_to_unicode(replace(expr, r"\|\-?\d+\|" => ""))
    return Formula{valtype(composition)}(
        expr, phreeqc_expr, unicode_expr, colored_formula(unicode_expr), composition, charge
    )
end

"""
    Formula(composition::AbstractDict{Symbol,T}, charge=0; order=ATOMIC_ORDER) where {T<:Number} -> Formula{T}

Construct a `Formula` from an explicit composition mapping.

# Arguments

  - `composition`: mapping of Symbol to numeric coefficient.
  - `charge`: explicit integer charge (default 0).
  - `order`: atomic ordering used for serialization.

Charge placeholder keys (`:Zz`, `:Zz⁺`, `:e`, `:e⁻`) are removed from the stored composition
and used to compute the formal charge when `charge == 0`.

# Examples

```jldoctest
julia> f = Formula(OrderedDict(:Ca=>1, :C=>1, :O=>3));

julia> expr(f)
"CaCO3"
```
"""
function Formula(
    composition::AbstractDict{Symbol,T}, charge=0; order=ATOMIC_ORDER
) where {T<:Number}
    charge_symbols = [:Zz, :Zz⁺, :e, :e⁻]
    filtered_keys = setdiff(keys(composition), charge_symbols)
    sorted_keys = sort(collect(filtered_keys); by=k -> findfirst(==(k), order))

    expr_parts = String[]
    uni_parts = String[]
    col_expr_parts = String[]
    for k in sorted_keys
        v = stoich_coef_round(composition[k])
        if !iszero(v)
            strv0 = get(dict_frac_unicode, v, string(v))
            strv = replace(strv0, " " => "", "*" => "")
            strvuni = strv
            colstrv = strv
            if occursin("+", strv0) || occursin("-", strv0) || occursin("*", strv0)
                strv = "(" * strv * ")"
            end
            if any(x -> x ∉ keys(dict_all_normal_to_sub), strv)
                strvuni = strv
                colstrv = string(COL_STOICH_INT(strv))
            else
                strvuni = all_normal_to_sub(strvuni)
                colstrv = string(COL_STOICH_INT(strvuni))
            end
            push!(expr_parts, string(k) * (isone(v) ? "" : strv))
            push!(uni_parts, string(k) * (isone(v) ? "" : strvuni))
            push!(col_expr_parts, string(k) * (isone(v) ? "" : colstrv))
        end
    end
    expr = join(expr_parts, "")
    uni = join(uni_parts, "")
    col_expr = join(col_expr_parts, "")

    if iszero(charge)
        charge = get(composition, :Zz, 0) + get(composition, :Zz⁺, 0)
        charge -= get(composition, :e, 0) + get(composition, :e⁻, 0)
    end
    if !iszero(charge)
        sign = charge < 0 ? "-" : "+"
        abscharge = abs(charge)
        strch = isone(abscharge) ? "" : string(abscharge)
        expr *= sign * strch
        uni *= sign * strch
        col_expr *= string(COL_CHARGE(normal_to_super(sign * strch)))
    end

    newcomposition = OrderedDict(
        k => v for (k, v) in composition if k ∉ charge_symbols && !iszero(v)
    )

    return Formula{T}(expr, unicode_to_phreeqc(expr), uni, col_expr, newcomposition, charge)
end

"""
    Formula(f::Formula) -> Formula

Copy constructor: return a new `Formula` built from `f`'s composition.
"""
function Formula(f::Formula)
    return Formula(composition(f))
end

"""
    Base.getindex(f::Formula{T}, i::Symbol) where {T} -> T

Return the stoichiometric coefficient associated with symbol `i`.
If `i` is not present, return `zero(T)` where `T` is the formula's coefficient type.

# Examples

```jldoctest
julia> f = Formula("H2O");

julia> f[:H]
2

julia> f[:C]
0
```
"""
function Base.getindex(f::Formula{T}, i::Symbol) where {T}
    coef = get(composition(f), i, nothing)
    if isnothing(coef)
        return zero(T)
    end
    return coef
end

"""
    Base.length(f::Formula) -> Int

Return the number of distinct element symbols in the formula composition.

# Examples

```jldoctest
julia> length(Formula("(CaO)1.25(SiO2)1(Al2O3)0.125(Na2O)0.25(H2O)1.375"))
6
```
"""
Base.length(f::Formula) = length(composition(f))

"""
    Base.isequal(f1::Formula, f2::Formula) -> Bool

Two formulas are equal if their compositions and formal charges are equal.

# Examples

```jldoctest
julia> Formula("H2O") == Formula(OrderedDict(:H=>2, :O=>1))
true
```
"""
function Base.isequal(f1::Formula, f2::Formula)
    isequal(composition(f1), composition(f2)) && isequal(charge(f1), charge(f2))
end
==(f1::Formula, f2::Formula) = isequal(f1, f2)

"""
    Base.hash(f::Formula, h::UInt) -> UInt

Hash a `Formula` using its composition and charge for stable use in collections.
"""
Base.hash(f::Formula, h::UInt) = hash(composition(f), hash(charge(f), h))

"""
    Base.show(io::IO, f::Formula)

Concise single-line representation for `Formula` objects, joining available textual forms.
"""
function Base.show(io::IO, f::Formula)
    print(io, join(unique!([expr(f), phreeqc(f), unicode(f)]), " ◆ "))
end

"""
    print_formula(io::IO, f::Formula, title::String, pad::Int)

Helper to print a titled, padded multi-field representation of a Formula.
Used by the MIME text/plain `show` method.

# Arguments

  - `io`: I/O stream.
  - `f`: formula to print.
  - `title`: section title.
  - `pad`: left-padding width.
"""
function print_formula(io::IO, f::Formula, title::String, pad::Int)
    println(
        io,
        lpad(title, pad),
        ": ",
        join(unique!([expr(f), phreeqc(f), unicode(f), colored(f)]), " ◆ "),
    )
end

"""
    Base.show(io::IO, ::MIME"text/plain", f::Formula)

Detailed multi-line pretty-printing used by the REPL.
Shows type, formula, composition and charge.
"""
function Base.show(io::IO, ::MIME"text/plain", f::Formula)
    pad = 11
    println(io, typeof(f))
    print_formula(io, f, "formula", pad)
    println(
        io,
        lpad("composition", pad),
        ": ",
        join(["$k => $v" for (k, v) in composition(f)], ", "),
    )
    println(io, lpad("charge", pad), ": ", charge(f))
end

"""
    *(f::Formula, x::T) where {T<:Number} -> Formula

Multiply all stoichiometric coefficients of `f` by scalar `x`.
"""
function *(f::Formula, x::T) where {T<:Number}
    composition = OrderedDict(k => x * v for (k, v) in f.composition)
    return Formula(composition, f.charge)
end

"""
    /(f::Formula, x::T) where {T<:Number} -> Formula

Divide all stoichiometric coefficients of `f` by scalar `x`.
"""
function /(f::Formula, x::T) where {T<:Number}
    composition = OrderedDict(k => v / x for (k, v) in f.composition)
    return Formula(composition, f.charge)
end

"""
    //(f::Formula, x::T) where {T<:Number} -> Formula

Produce rational coefficients by dividing `f`'s coefficients by `x` (rational result).
"""
function //(f::Formula, x::T) where {T<:Number}
    composition = OrderedDict(k => v // x for (k, v) in f.composition)
    return Formula(composition, f.charge)
end

"""
    +(f::Formula, atom::AtomGroup) -> Formula

Add an `AtomGroup` to a Formula (adjust stoichiometric coefficient for the atom).
"""
function +(f::Formula, atom::AtomGroup)
    composition = copy(f.composition)
    composition[atom.sym] = get(f.composition, atom.sym, 0) + atom.coef
    return Formula(composition, f.charge)
end

"""
    +(a::AtomGroup{T}, b::AtomGroup{S}) where {T,S} -> Formula

Combine two `AtomGroup` values into a `Formula`. If the symbols are equal the
result is a singleton composition with summed coefficients, otherwise both are
included.

# Examples

```jldoctest
julia> result = AtomGroup(:H, 2) + AtomGroup(:H, 1);

julia> composition(result)[:H]
3
```
"""
function +(a::AtomGroup{T}, b::AtomGroup{S}) where {T,S}
    composition = if a.sym == b.sym
        OrderedDict{Symbol,promote_type(T, S)}(a.sym => a.coef + b.coef)
    else
        Dict{Symbol,promote_type(T, S)}(a.sym => a.coef, b.sym => b.coef)
    end
    return Formula(composition, 0)
end

"""
    +(a::AtomGroup, b::Symbol) -> Formula

Convenience: add an AtomGroup and a Symbol (converted to an AtomGroup of coef 1).
"""
+(a::AtomGroup, b::Symbol) = a + AtomGroup(b)

"""
    Base.convert(T::Type{<:Number}, f::Formula) -> Formula{T}

Convert the stoichiometric coefficient type of `f` to numeric type `T`.
"""
function Base.convert(T::Type{<:Number}, f::Formula)
    newcomposition = OrderedDict(k => convert(T, v) for (k, v) in composition(f))
    return Formula{T}(
        expr(f), phreeqc(f), unicode(f), colored(f), newcomposition, charge(f)
    )
end

"""
    apply(func::Function, f::Formula, args...; kwargs...) -> Formula

Element-wise apply `func` to all numeric components of `f` and to its charge.
Quantities are handled, attempting to preserve dimensions when possible.

# Examples

```jldoctest
julia> result = apply(x -> x*2, Formula("H2O"));

julia> result[:H]
4
```
"""
function apply(func::Function, f::Formula, args...; kwargs...)
    tryfunc(v) =
        if v isa Quantity
            (
                try
                    func(ustrip(v), args...; kwargs...) *
                    func(dimension(v), args...; kwargs...)
                catch
                    try
                        func(ustrip(v), args...; kwargs...) * dimension(v)
                    catch
                        v
                    end
                end
            )
        else
            (
                try
                    func(v, args...; kwargs...)
                catch
                    v
                end
            )
        end
    newcomposition = OrderedDict(k => tryfunc(v) for (k, v) in composition(f))
    return Formula{valtype(newcomposition)}(
        expr(f), phreeqc(f), unicode(f), colored(f), newcomposition, tryfunc(charge(f))
    )
end

"""
    check_mendeleev(f::Formula) -> Bool

Validate that all element symbols in `f` exist in the package `elements`
registry. Returns `true` when valid; otherwise throws an informative error.

# Examples

```jldoctest
julia> check_mendeleev(Formula("NaCl"))
true
```
"""
function check_mendeleev(f::Formula)
    nonatoms = filter(x -> x ∉ keys(elements.bysymbol) && x != :Zz, keys(composition(f)))
    # isempty(nonatoms) ? true : error("Invalid atoms $(join(nonatoms," ")) in $f")
    return isempty(nonatoms)
end
