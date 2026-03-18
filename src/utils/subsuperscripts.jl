const dict_super_to_normal = OrderedDict{Char, Char}(
    '⁰' => '0',
    '¹' => '1',
    '²' => '2',
    '³' => '3',
    '⁴' => '4',
    '⁵' => '5',
    '⁶' => '6',
    '⁷' => '7',
    '⁸' => '8',
    '⁹' => '9',
    '⁺' => '+',
    '⁻' => '-',
    '.' => '.',
)

const dict_normal_to_super = OrderedDict{Char, Char}(
    '0' => '⁰',
    '1' => '¹',
    '2' => '²',
    '3' => '³',
    '4' => '⁴',
    '5' => '⁵',
    '6' => '⁶',
    '7' => '⁷',
    '8' => '⁸',
    '9' => '⁹',
    '+' => '⁺',
    '-' => '⁻',
    '.' => '.',
)

const dict_sub_to_normal = OrderedDict{Char, Char}(
    '₀' => '0',
    '₁' => '1',
    '₂' => '2',
    '₃' => '3',
    '₄' => '4',
    '₅' => '5',
    '₆' => '6',
    '₇' => '7',
    '₈' => '8',
    '₉' => '9',
    '.' => '.',
)

const dict_normal_to_sub = OrderedDict{Char, Char}(
    '0' => '₀',
    '1' => '₁',
    '2' => '₂',
    '3' => '₃',
    '4' => '₄',
    '5' => '₅',
    '6' => '₆',
    '7' => '₇',
    '8' => '₈',
    '9' => '₉',
    '+' => '₊',
    '-' => '₋',
    '.' => '.',
)

const dict_all_normal_to_sub = OrderedDict{Char, Char}(
    '0' => '₀',
    '1' => '₁',
    '2' => '₂',
    '3' => '₃',
    '4' => '₄',
    '5' => '₅',
    '6' => '₆',
    '7' => '₇',
    '8' => '₈',
    '9' => '₉',
    '+' => '₊',
    '-' => '₋',
    '=' => '₌',
    '(' => '₍',
    ')' => '₎',
    '.' => '.',
    'a' => 'ₐ',
    'e' => 'ₑ',
    'o' => 'ₒ',
    'x' => 'ₓ',
    'h' => 'ₕ',
    'k' => 'ₖ',
    'l' => 'ₗ',
    'm' => 'ₘ',
    'n' => 'ₙ',
    'p' => 'ₚ',
    's' => 'ₛ',
    't' => 'ₜ',
    '∂' => 'ₔ',
    'β' => 'ᵦ',
    'γ' => 'ᵧ',
    'ρ' => 'ᵨ',
    'φ' => 'ᵩ',
    'χ' => 'ᵪ',
    '*' => '*',
    '/' => '/',
)

const dict_frac_unicode = OrderedDict(
    1 // 4 => "¼",
    1 // 2 => "½",
    3 // 4 => "¾",
    1 // 7 => "⅐",
    1 // 9 => "⅑",
    1 // 10 => "⅒",
    1 // 3 => "⅓",
    2 // 3 => "⅔",
    1 // 5 => "⅕",
    2 // 5 => "⅖",
    3 // 5 => "⅗",
    4 // 5 => "⅘",
    1 // 6 => "⅙",
    5 // 6 => "⅚",
    1 // 8 => "⅛",
    3 // 8 => "⅜",
    5 // 8 => "⅝",
    7 // 8 => "⅞",
)

const dict_unicode_frac = OrderedDict(
    '¼' => 1 // 4,
    '½' => 1 // 2,
    '¾' => 3 // 4,
    '⅐' => 1 // 7,
    '⅑' => 1 // 9,
    '⅒' => 1 // 10,
    '⅓' => 1 // 3,
    '⅔' => 2 // 3,
    '⅕' => 1 // 5,
    '⅖' => 2 // 5,
    '⅗' => 3 // 5,
    '⅘' => 4 // 5,
    '⅙' => 1 // 6,
    '⅚' => 5 // 6,
    '⅛' => 1 // 8,
    '⅜' => 3 // 8,
    '⅝' => 5 // 8,
    '⅞' => 7 // 8,
)

"""
    issuperscript(c::Char) -> Bool

Return whether `c` is a numeric superscript or ⁺/⁻.

# Examples

```julia
julia> issuperscript('²')
true

julia> issuperscript('a')
false
```
"""
issuperscript(c::Char) = c in keys(dict_super_to_normal)

"""
    issubscript(c::Char) -> Bool

Return whether `c` is a numeric subscript.

# Examples

```julia
julia> issubscript('₂')
true

julia> issubscript('2')
false
```
"""
issubscript(c::Char) = c in keys(dict_sub_to_normal)

"""
    super_to_normal(s::AbstractString) -> String

Convert all numeric superscripts or ⁺/⁻ in `s` to normal line.

# Examples

```julia
julia> super_to_normal("Ca²⁺")
"Ca2+"
```
"""
super_to_normal(s::AbstractString) = replace(s, dict_super_to_normal...)

"""
    normal_to_super(s::AbstractString) -> String

Convert all normal characters or +/- in `s` to numeric superscripts.

# Examples

```julia
julia> normal_to_super("2+")
"²⁺"
```
"""
normal_to_super(s::AbstractString) = replace(s, dict_normal_to_super...)

"""
    sub_to_normal(s::AbstractString) -> String

Convert all numeric subscripts in `s` to normal line.

# Examples

```julia
julia> sub_to_normal("H₂O")
"H2O"
```
"""
sub_to_normal(s::AbstractString) = replace(s, dict_sub_to_normal...)

"""
    normal_to_sub(s::AbstractString) -> String

Convert all normal characters in `s` to numeric subscripts.

# Examples

```julia
julia> normal_to_sub("H2O")
"H₂O"
```
"""
normal_to_sub(s::AbstractString) = replace(s, dict_normal_to_sub...)

"""
    all_normal_to_sub(s::AbstractString) -> String

Convert all normal characters (including letters and operators) to subscripts.
"""
all_normal_to_sub(s::AbstractString) = replace(s, dict_all_normal_to_sub...)

"""
    subscriptnumber(i::Integer) -> String

Convert an integer to its Unicode subscript representation.

# Examples

```julia
julia> subscriptnumber(42)
"₄₂"

julia> subscriptnumber(-3)
"₋₃"
```
"""
function subscriptnumber(i::Integer)
    if i < 0
        c = [Char(0x208B)]
    else
        c = []
    end
    for d in reverse(digits(abs(i)))
        push!(c, Char(0x2080 + d))
    end
    return join(c)
end

"""
    superscriptnumber(i::Integer) -> String

Convert an integer to its Unicode superscript representation.

# Examples

```julia
julia> superscriptnumber(42)
"⁴²"

julia> superscriptnumber(-2)
"⁻²"
```
"""
function superscriptnumber(i::Integer)
    if i < 0
        c = [Char(0x207B)]
    else
        c = []
    end
    for d in reverse(digits(abs(i)))
        if d == 0
            push!(c, Char(0x2070))
        end
        if d == 1
            push!(c, Char(0x00B9))
        end
        if d == 2
            push!(c, Char(0x00B2))
        end
        if d == 3
            push!(c, Char(0x00B3))
        end
        if d > 3
            push!(c, Char(0x2070 + d))
        end
    end
    return join(c)
end

"""
    from_subscriptnumber(s::String) -> Int

Parse a Unicode subscript number string to an integer.

# Examples

```julia
julia> from_subscriptnumber("₄₂")
42

julia> from_subscriptnumber("₋₃")
-3
```
"""
function from_subscriptnumber(s::String)
    chars = collect(s)
    negative = !isempty(chars) && chars[1] == Char(0x208B)
    if negative
        chars = chars[2:end]
    end
    value = 0
    for c in chars
        digit = Int(c) - 0x2080
        value = value * 10 + digit
    end
    return negative ? -value : value
end

"""
    from_superscriptnumber(s::String) -> Int

Parse a Unicode superscript number string to an integer.

# Examples

```julia
julia> from_superscriptnumber("⁴²")
42

julia> from_superscriptnumber("⁻²")
-2
```
"""
function from_superscriptnumber(s::String)
    chars = collect(s)
    negative = !isempty(chars) && chars[1] == Char(0x207B)
    if negative
        chars = chars[2:end]
    end
    value = 0
    for c in chars
        if c == Char(0x2070)
            digit = 0
        elseif c == Char(0x00B9)
            digit = 1
        elseif c == Char(0x00B2)
            digit = 2
        elseif c == Char(0x00B3)
            digit = 3
        else
            digit = Int(c) - 0x2070
        end
        value = value * 10 + digit
    end
    return negative ? -value : value
end
