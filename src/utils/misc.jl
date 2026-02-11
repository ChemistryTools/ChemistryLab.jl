const COL_CHARGE = crayon"cyan bold"
const COL_PAR = crayon"magenta bold"
const COL_STOICH_INT = crayon"red bold"
const COL_STOICH_EXT = crayon"yellow bold"

"""
    print_title(title; crayon=:none, indent="", style=:none)

Print a formatted title with optional styling and indentation.

# Arguments

  - `title`: title string to display.
  - `crayon`: Crayon for colored output (default `:none`).
  - `indent`: left indentation string (default `""`).
  - `style`: formatting style - `:none`, `:underline`, or `:box` (default `:none`).

# Examples

```julia
julia> print_title("Test Title"; style=:none)
Test Title

julia> print_title("Section"; style=:underline, indent="  ")
  Section
  ───────
```
"""
function print_title(title; crayon=:none, indent="", style=:none)
    draw = crayon != :none ? x -> println(crayon(x)) : println
    if style == :underline
        width = length(title)
        draw(indent * "$title")
        draw(indent * "─"^width)
    elseif style == :box
        width = length(title) + 4
        draw(indent * "┌" * "─"^width * "┐")
        draw(indent * "│  $title  │")
        draw(indent * "└" * "─"^width * "┘")
    else
        draw(indent * "$title")
    end
    print(crayon"reset")
end

"""
    root_type(T) -> Type

Get the root type wrapper from a potentially parameterized type.
"""
root_type(T) = T isa UnionAll ? root_type(T.body) : T.name.wrapper

"""
    safe_ustrip(unit::UnionAbstractQuantity, q::UnionAbstractQuantity) -> Number
    safe_ustrip(unit::UnionAbstractQuantity, q) -> Number

Safely remove units from `q`. If `q` is a quantity the returned numeric value
is expressed in `unit`. If `q` is a plain numeric value, it is returned
unchanged.

# Arguments

  - `unit` : target unit to express `q` in (ignored if `q` has no units).
  - `q` : quantity or numeric value.

# Returns

  - Numeric value of `q` expressed in `unit` (when `q` had units) or `q`
    unchanged (when `q` had no units).

# Examples

```julia
julia> safe_ustrip(u"m", 5u"m")
5.0

julia> safe_ustrip(u"cm", 1u"m")
100.0

julia> safe_ustrip(u"m", 3.0)
3.0
```
"""
safe_ustrip(unit::UnionAbstractQuantity, q::UnionAbstractQuantity) = ustrip(unit, q)

safe_ustrip(::UnionAbstractQuantity, q) = ustrip(q)

"""
    safe_uconvert(qout::UnionAbstractQuantity{<:Any, <:AbstractSymbolicDimensions}, q::UnionAbstractQuantity{<:Any, <:Dimensions}) -> UnionAbstractQuantity
    safe_uconvert(qout::UnionAbstractQuantity{<:Any, <:AbstractSymbolicDimensions}, q) -> Any

Convert `q` to the unit type represented by `qout` using `uconvert` when `q`
is a quantity with compatible dimensions. If `q` is a plain numeric value,
return `q` unchanged.

# Arguments

  - `qout` : unitful quantity target (e.g., u"m", u"kg").
  - `q` : a quantity or a plain numeric value.

# Returns

  - Converted quantity (of type `qout`) or the original `q` when `q` has no units.

# Examples

```julia
julia> safe_uconvert(us"m", 5u"cm")
0.05 m

julia> safe_uconvert(us"m", 3.0)
3.0
```
"""
safe_uconvert(qout::UnionAbstractQuantity{<:Any, <:AbstractSymbolicDimensions}, q::UnionAbstractQuantity{<:Any, <:Dimensions}) = uconvert(qout, q)

safe_uconvert(::UnionAbstractQuantity{<:Any, <:AbstractSymbolicDimensions}, q) = q

safe_uparse(x::AbstractString) = uparse(x)

safe_uparse(x::AbstractQuantity) = x

# force_uconvert(qout::UnionAbstractQuantity{<:Any, <:AbstractSymbolicDimensions}, q::UnionAbstractQuantity{<:Any, <:Dimensions}) = safe_uconvert(qout, q)

force_uconvert(qout::UnionAbstractQuantity, q) = safe_ustrip(qout, q)*qout
