*(x::Num, y::AbstractQuantity) = Quantity(x * y.value, y.dimensions)
*(x::AbstractQuantity, y::Num) = *(y, x)

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

const ATOMIC_ORDER = [
    :Ca,
    :Na,
    :K,
    :Mg,
    :Sr,
    :Ba,
    :Al,
    :Fe,
    :Ti,
    :Mn,
    :Cr,
    :Si,
    :C,
    :H,
    :N,
    :S,
    :O,
    :P,
    :B,
    :F,
    :Cl,
    :Br,
    :I,
    :Zz,
]

const cement_to_mendeleev = [
    :C => OrderedDict(:Ca => 1, :O => 1),
    :M => OrderedDict(:Mg => 1, :O => 1),
    :S => OrderedDict(:Si => 1, :O => 2),
    :A => OrderedDict(:Al => 2, :O => 3),
    :F => OrderedDict(:Fe => 2, :O => 3),
    :K => OrderedDict(:K => 2, :O => 1),
    :N => OrderedDict(:Na => 2, :O => 1),
    :P => OrderedDict(:P => 2, :O => 5),
    :T => OrderedDict(:Ti => 1, :O => 2),
    :C̄ => OrderedDict(:C => 1, :O => 2),
    :S̄ => OrderedDict(:S => 1, :O => 3),
    :N̄ => OrderedDict(:N => 1, :O => 3),
    :H => OrderedDict(:H => 2, :O => 1),
]

const OXIDE_ORDER = collect(first.(cement_to_mendeleev))

const dict_super_to_normal = OrderedDict{Char,Char}(
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

const dict_normal_to_super = OrderedDict{Char,Char}(
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

const dict_sub_to_normal = OrderedDict{Char,Char}(
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

const dict_normal_to_sub = OrderedDict{Char,Char}(
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

const dict_all_normal_to_sub = OrderedDict{Char,Char}(
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
    fwd_arrows

Collection of forward arrow symbols used in chemical reaction notation.

Contains symbols representing forward reaction directions, including:

  - Simple arrows: '>', '→'
  - Various Unicode arrows with different styles and weights
  - Specialized arrows for reaction notation

Used to define reaction directionality from reactants to products.
"""
const fwd_arrows = ['>', '→', '↣', '↦', '⇾', '⟶', '⟼', '⥟', '⥟', '⇀', '⇁', '⇒', '⟾']

"""
    bwd_arrows

Collection of backward arrow symbols used in chemical reaction notation.

Contains symbols representing reverse reaction directions, including:

  - Simple arrows: '<', '←'
  - Various Unicode arrows with different styles and weights
  - Specialized arrows for reverse reaction notation

Used to define reaction directionality from products to reactants.
"""
const bwd_arrows = ['<', '←', '↢', '↤', '⇽', '⟵', '⟻', '⥚', '⥞', '↼', '↽', '⇐', '⟽']

"""
    double_arrows

Collection of double arrow symbols representing equilibrium reactions.

Contains symbols representing:

  - Simple equilibrium: '↔'
  - Various Unicode double arrows
  - Specialized equilibrium symbols
  - Bidirectional reaction indicators

Used to denote reversible reactions and equilibrium states.
"""
const double_arrows = ['↔', '⟷', '⇄', '⇆', '⇌', '⇋', '⇔', '⟺']

"""
    pure_rate_arrows

Collection of specialized arrows for rate-based reaction notation.

Contains symbols commonly used to represent:

  - Reaction rates
  - Kinetic directions
  - Specialized reaction mechanisms

These are often used in more advanced chemical kinetics notation.
"""
const pure_rate_arrows = ['⇐', '⟽', '⇒', '⟾', '⇔', '⟺']

"""
    equal_signs

Collection of equality signs used in chemical reaction equations.

Contains various forms of equality operators including:

  - Standard equals sign: '='
  - Definition operators: '≔'
  - Specialized equality symbols
  - Assignment operators

Used to separate reactants from products in balanced equations.
"""
const equal_signs = ['=', '≔', '⩴', '≕']

const EQUAL_REACTION = vcat(
    fwd_arrows, bwd_arrows, double_arrows, pure_rate_arrows, equal_signs
)
const EQUAL_REACTION_SET = Set(EQUAL_REACTION)

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

"""
    parse_symbol_num(s::AbstractString) -> NamedTuple

Parse a symbol with optional numeric index into components.

# Arguments

  - `s`: string potentially containing a base symbol and numeric index.

# Returns

  - NamedTuple with `base` (symbol part), `index` (numeric part, or -100 if none),
    and `convert_func` (function to convert index back to original format).

# Examples

```julia
julia> parse_symbol_num("T₂₉₈")
(base = "T", index = 298, convert_func = subscriptnumber)

julia> parse_symbol_num("Ca")
(base = "Ca", index = -100, convert_func = identity)
```
"""
function parse_symbol_num(s::AbstractString)
    chars = collect(s)
    n = length(chars)
    pos = n + 1
    for i in 1:n
        c = chars[i]
        if '0' <= c <= '9'
            pos = i
            break
        end
        if 0x2080 <= Int(c) <= 0x2089
            pos = i
            break
        end
        if c == Char(0x208B)  # signe moins subscript
            pos = i
            break
        end
        if c in (Char(0x2070), Char(0x00B9), Char(0x00B2), Char(0x00B3)) ||
            (0x2074 <= Int(c) <= 0x2079)
            pos = i
            break
        end
    end
    if pos == n + 1
        return (base=s, index=-100, convert_func=identity)
    end
    radical = join(chars[1:(pos - 1)])
    number_chars = chars[pos:end]
    firstnum = number_chars[1]
    number = join(number_chars)
    if '0' <= firstnum <= '9'
        nature = identity
        number = parse(Int, number)
    elseif 0x2080 <= Int(firstnum) <= 0x2089 || firstnum == Char(0x208B)
        nature = subscriptnumber
        number = from_subscriptnumber(number)
    elseif firstnum in (Char(0x2070), Char(0x00B9), Char(0x00B2), Char(0x00B3)) ||
        (0x2074 <= Int(firstnum) <= 0x2079)
        nature = superscriptnumber
        number = from_superscriptnumber(number)
    else
        nature = identity
        number = 0
    end
    return (base=radical, index=number, convert_func=nature)
end

"""
    extract_vars(expr) -> Vector{Symbol}

Extract all variable symbols from a symbolic expression.

# Arguments

  - `expr`: Julia expression (Expr or Symbol).

# Returns

  - Vector of unique Symbol variables found in the expression.
"""
function extract_vars(expr)
    vars = Set{Symbol}()
    function _extract(e, is_func=false)
        if isa(e, Symbol)
            if !is_func
                push!(vars, e)
            end
        elseif isa(e, Expr)
            if e.head == :call
                _extract(e.args[1], true)
                for arg in e.args[2:end]
                    _extract(arg, false)
                end
            else
                for arg in e.args
                    _extract(arg, false)
                end
            end
        end
    end
    _extract(expr, false)
    return collect(vars)
end

"""
    root_type(T) -> Type

Get the root type wrapper from a potentially parameterized type.
"""
root_type(T) = T isa UnionAll ? root_type(T.body) : T.name.wrapper

"""
    stoich_coef_round(x::T; tol=1e-4) where {T<:Real} -> Union{Int, Rational, Float64}
    stoich_coef_round(x) -> Any

Round stoichiometric coefficients to integer, rational, or float representation.

# Arguments

  - `x`: numeric value to round.
  - `tol`: tolerance for rounding decisions (default 1e-4).

# Returns

  - Integer if close to a whole number, Rational if a simple fraction (denominator < 10),
    or Float64 rounded to 5 digits otherwise. Non-numeric inputs are returned unchanged.

# Examples

```julia
julia> stoich_coef_round(2.0001)
2

julia> stoich_coef_round(0.3333)
1//3

julia> stoich_coef_round(3.14159)
3.14159
```
"""
function stoich_coef_round(x::T; tol=1e-4) where {T<:Real}
    try
        if isapprox(x, round(x); atol=tol)
            return Int(round(x))
        end

        rat = rationalize(x; tol=tol)
        if isapprox(x, float(rat); atol=tol)
            if 1 < denominator(rat) < 10
                return rat
            end
        end

        return round(x; digits=5)
    catch e
        return x
    end
end

stoich_coef_round(x) = x

"""
    phreeqc_to_unicode(s::AbstractString) -> String

Convert a PHREEQC formula string to Unicode representation with subscripts and superscripts.

# Arguments

  - `s`: PHREEQC-formatted chemical formula string.

# Returns

  - Unicode-formatted string with subscript coefficients and superscript charges.

# Examples

```julia
julia> phreeqc_to_unicode("Ca+2")
"Ca²⁺"

julia> phreeqc_to_unicode("SO4-2")
"SO₄²⁻"
```
"""
function phreeqc_to_unicode(s::AbstractString)
    chars = collect(s)

    ind_sign = findall(
        i -> (chars[i] == '+' || chars[i] == '-') && (i > 1 && chars[i - 1] != ' '),
        1:length(chars),
    )

    for i in ind_sign
        sign = chars[i]
        j = i
        while j < length(chars) && (isnumeric(chars[j + 1]) || chars[j + 1] == '.')
            chars[j] = dict_normal_to_super[chars[j + 1]]
            j += 1
        end
        chars[j] = dict_normal_to_super[sign]
    end

    s = join(chars)

    s = replace(s, r"-?\d+\.?\d*" => x -> string(stoich_coef_round(parse(Float64, x))))

    matches = collect(eachmatch(r"(\d+)\/\/(\d+)", s))
    for m in reverse(matches)
        num = tryparse(Int, m.captures[1])
        den = tryparse(Int, m.captures[2])
        rat = stoich_coef_round(num // den)
        replacement = get(dict_frac_unicode, rat, string(rat))
        start_idx = m.offset
        end_idx = start_idx + sizeof(m.match) - 1
        prefix = if start_idx > first(eachindex(s))
            s[first(eachindex(s)):prevind(s, start_idx)]
        else
            ""
        end
        suffix = end_idx < sizeof(s) ? s[(end_idx + 1):end] : ""
        s = prefix * replacement * suffix
    end

    chars = collect(s)

    ind_sign = findall(
        i ->
            chars[i] in keys(dict_normal_to_sub) &&
            i > 1 &&
            chars[i - 1] != ' ' &&
            !(chars[i - 1] in keys(dict_normal_to_sub)),
        1:length(chars),
    )

    for i in ind_sign
        j = i
        while j <= length(chars) && chars[j] in keys(dict_normal_to_sub)
            chars[j] = dict_normal_to_sub[chars[j]]
            j += 1
        end
    end

    return join(chars)
end

"""
    merge_upper_lower(graphemes::Vector{<:AbstractString}) -> Vector{String}

Merge consecutive graphemes where an uppercase letter is followed by a lowercase letter.

# Arguments

  - `graphemes`: vector of grapheme strings.

# Returns

  - Vector with merged element symbols (e.g., ["C", "a"] becomes ["Ca"]).
"""
function merge_upper_lower(graphemes::Vector{<:AbstractString})
    result = String[]
    i = 1
    while i <= length(graphemes)
        current = graphemes[i]
        if i < length(graphemes)
            last_char = current[end]
            next_first_char = graphemes[i + 1][1]
            if isuppercase(last_char) && islowercase(next_first_char)
                current *= graphemes[i + 1]
                i += 1
            end
        end
        push!(result, current)
        i += 1
    end
    return result
end

"""
    unicode_to_phreeqc(s::AbstractString) -> String

Convert a Unicode formula string back to PHREEQC format.

# Arguments

  - `s`: Unicode-formatted chemical formula.

# Returns

  - PHREEQC-formatted string with plain text charges and coefficients.

# Examples

```julia
julia> unicode_to_phreeqc("Ca²⁺")
"Ca+2"

julia> unicode_to_phreeqc("SO₄²⁻")
"SO4-2"
```
"""
function unicode_to_phreeqc(s::AbstractString)
    chars = collect(s)

    ind_sign = findall(k -> k == '⁺' || k == '⁻', chars)
    for i in ind_sign
        sign = chars[i]
        j = i
        while j > 1 && chars[j - 1] in keys(dict_super_to_normal)
            chars[j] = chars[j - 1]
            j -= 1
        end
        chars[j] = sign
    end
    s = super_to_normal(join(chars))

    s = replace(s, r"[₀₁₂₃₄₅₆₇₈₉]" => c -> dict_sub_to_normal[c[1]])

    pattern = r"(\d+(\.\d+)?|)([¼½¾⅐⅑⅒⅓⅔⅕⅖⅗⅘⅙⅚⅛⅜⅝⅞])"
    matches = collect(eachmatch(pattern, s))
    for m in reverse(matches)
        float_part = isempty(m.captures[1]) ? 0.0 : parse(Float64, m.captures[1])
        sum_value = float_part + dict_unicode_frac[m.captures[3][1]]
        replacement = string(stoich_coef_round(sum_value))
        start_idx = m.offset
        end_idx = start_idx + sizeof(m.match) - 1
        prefix = if start_idx > first(eachindex(s))
            s[first(eachindex(s)):prevind(s, start_idx)]
        else
            ""
        end
        suffix = end_idx < sizeof(s) ? s[(end_idx + 1):end] : ""
        s = prefix * replacement * suffix
    end

    return s
end

"""
    colored_formula(s::AbstractString; colorcharge=true) -> String

Generate a terminal-colored representation of a chemical formula.

# Arguments

  - `s`: formula string.
  - `colorcharge`: if true, color the charge portion (default true).

# Returns

  - String with ANSI color codes for terminal display.
"""
function colored_formula(s::AbstractString; colorcharge=true)
    superscript_digits = ["⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹", "⁺", "⁻"]

    colored_graph = merge_upper_lower(collect(graphemes(s)))
    ind_sign = findlast(k -> k == "⁺" || k == "⁻" || k == "+" || k == "-", colored_graph)
    if isnothing(ind_sign)
        ind_sign = length(colored_graph) + 1
    end

    idx_sign = Int[]
    idx_atoms = Int[]
    idx_par = Int[]
    for (i, c) in enumerate(colored_graph)
        if colorcharge && (i >= ind_sign || c in superscript_digits)
            push!(idx_sign, i)
        end
        if Symbol(c) in ATOMIC_ORDER || Symbol(c) in OXIDE_ORDER
            push!(idx_atoms, i)
        end
        if c in ["(", ")", "[", "]", "{", "}", "@", "|"]
            push!(idx_par, i)
        end
    end
    idx_stoich = setdiff(1:length(colored_graph), union(idx_sign, idx_atoms, idx_par))
    colored_graph[idx_sign] .= string.(COL_CHARGE.(colored_graph[idx_sign]))
    colored_graph[idx_par] .= string.(COL_PAR.(colored_graph[idx_par]))
    colored_graph[idx_stoich] .= string.(COL_STOICH_INT.(colored_graph[idx_stoich]))

    return join(colored_graph)
end

"""
    parse_formula(formula::AbstractString) -> OrderedDict{Symbol,Number}

Parse a chemical formula string into an atomic composition dictionary.

# Arguments

  - `formula`: formula string (supports parentheses, brackets, rational/decimal coefficients).

# Returns

  - OrderedDict mapping element symbols to stoichiometric coefficients.

# Examples

```julia
julia> parse_formula("H2O")
OrderedDict{Symbol, Number} with 2 entries:
  :H => 2
  :O => 1

julia> parse_formula("Ca(OH)2")
OrderedDict{Symbol, Number} with 3 entries:
  :Ca => 1
  :O  => 2
  :H  => 2
```
"""
function parse_formula(formula::AbstractString)
    function safe_nextind(s::AbstractString, i::Integer, n::Integer=1)
        last_i = lastindex(s)
        idx = i
        for _ in 1:n
            if idx > last_i
                return last_i + 1
            end
            idx = nextind(s, idx)
        end
        return idx
    end

    formula = replace(formula, ":" => "", "{" => "(", "}" => ")", "[" => "(", "]" => ")")
    formula = replace(formula, r"\|\-?\d+\|" => "")
    formula = replace(formula, r"\|" => "")
    formula = unicode_to_phreeqc(String(formula))

    counts = OrderedDict{Symbol,Number}()

    i = firstindex(formula)
    while i <= lastindex(formula)
        c = formula[i]
        if c == '('
            depth = 1
            j = safe_nextind(formula, i)
            while j <= lastindex(formula) && depth > 0
                if formula[j] == '('
                    depth += 1
                elseif formula[j] == ')'
                    depth -= 1
                end
                j = safe_nextind(formula, j)
            end
            inner = formula[safe_nextind(formula, i):prevind(formula, j)]

            rest = j <= lastindex(formula) ? formula[j:end] : ""

            m = match(r"^([0-9]+//[0-9]+|[0-9]+(?:\.[0-9]+)?)", rest)
            factor = if (m === nothing)
                1
            else
                begin
                    s = m.match
                    factor =
                        occursin("//", s) ? parse(Rational{Int}, s) : parse(Float64, s)
                end
            end
            offset = (m === nothing) ? 0 : length(m.match)

            for (el, n) in parse_formula(inner)
                counts[el] = get(counts, el, 0) + n * factor
            end

            i = safe_nextind(formula, j, offset)
        else
            m = match(
                r"^(\p{Lu}[\p{Ll}\u0300-\u036F]?)(([0-9]+//[0-9]+)|([0-9]+(?:\.[0-9]+)?))?",
                formula[i:end],
            )
            if m !== nothing
                el, countstr = m.captures
                el = Symbol(el)
                cnt = if countstr === nothing || isempty(countstr)
                    1
                elseif occursin("//", countstr)
                    parse(Rational{Int}, countstr)
                else
                    stoich_coef_round(parse(Float64, countstr))
                end

                if cnt isa Rational && denominator(cnt) == 1
                    cnt = Int(numerator(cnt))
                elseif cnt isa Float64 && isinteger(cnt)
                    cnt = Int(cnt)
                end

                counts[el] = get(counts, el, 0) + cnt

                i = safe_nextind(formula, i, length(m.match))
            else
                i = safe_nextind(formula, i)
            end
        end
    end

    # T = promote_type(typeof.(stoich_coef_round.(values(counts)))...)

    return OrderedDict(k => stoich_coef_round(v) for (k, v) in counts)
end

"""
    extract_charge(formula::AbstractString) -> Int

Extract the formal charge from a chemical formula string.

# Arguments

  - `formula`: formula string with optional charge notation (e.g., "+2", "-", "3+").

# Returns

  - Integer charge value (0 if no charge present).

# Examples

```julia
julia> extract_charge("Ca+2")
2

julia> extract_charge("SO4-2")
-2

julia> extract_charge("H2O")
0
```
"""
function extract_charge(formula::AbstractString)
    m = match(r"([+-])([0-9]*)$", unicode_to_phreeqc(formula))
    if m === nothing
        return 0
    else
        sign = m.captures[1] == "+" ? 1 : -1
        val = m.captures[2] == "" ? 1 : parse(Int, m.captures[2])
        return sign * val
    end
end

"""
    calculate_molar_mass(atoms::AbstractDict{Symbol,T}) where {T<:Number} -> Quantity

Calculate the molar mass from an atomic composition dictionary.

# Arguments

  - `atoms`: dictionary mapping element symbols to stoichiometric coefficients.

# Returns

  - Molar mass as a Quantity in g/mol units.

# Examples

```julia
julia> calculate_molar_mass(OrderedDict(:H => 2, :O => 1))
18.01528 g mol⁻¹
```    # return sum(cnt * ustrip(elements[element].atomic_mass) for (element, cnt) in atoms if haskey(elements, element); init=0) * u"g/mol"
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

"""
    replace_graphemes(s::AbstractString, old_new::Pair...) -> String

Replace grapheme-level substrings in a Unicode string.

# Arguments

  - `s`: input string.
  - `old_new`: pairs of old => new grapheme replacements.

# Returns

  - String with replacements applied at grapheme level.
"""
function replace_graphemes(s::AbstractString, old_new::Pair...)
    gs = collect(graphemes(s))

    mapping = OrderedDict{String,String}()
    for pair in old_new
        mapping[string(pair.first)] = string(pair.second)
    end

    for i in eachindex(gs)
        if haskey(mapping, gs[i])
            gs[i] = mapping[gs[i]]
        end
    end
    return join(gs)
end

"""
    to_mendeleev(oxides::AbstractDict{Symbol,T}) where {T<:Number} -> OrderedDict{Symbol,Number}

Convert cement oxide notation to Mendeleev element composition.

# Arguments

  - `oxides`: dictionary mapping oxide symbols (C, S, A, etc.) to coefficients.

# Returns

  - OrderedDict mapping element symbols to stoichiometric coefficients.

# Examples

```julia
julia> to_mendeleev(OrderedDict(:C => 1, :S => 2))
OrderedDict{Symbol, Number} with 3 entries:
  :Ca => 1
  :O  => 5
  :Si => 2
```
"""
function to_mendeleev(oxides::AbstractDict{Symbol,T}) where {T<:Number}
    result = OrderedDict{Symbol,Number}()
    for (ox, coef) in oxides
        if ox ∉ [:Zz, :Zz⁺, :e, :e⁻]
            idx = findfirst(p -> p.first == ox, cement_to_mendeleev)
            idx !== nothing || error("$(ox) is not a valid oxide identifier")
            mend = cement_to_mendeleev[idx].second
            for (k, v) in mend
                result[k] = get(result, k, 0) + v * coef
            end
        end
    end
    return if length(result) > 0
        OrderedDict(k => stoich_coef_round(v) for (k, v) in result)
    else
        result
    end
end

"""
    parse_equation(equation::AbstractString) -> Tuple{OrderedDict{String,Real}, OrderedDict{String,Real}, Char}

Parse a chemical equation string into reactants, products, and equality sign.

# Arguments

  - `equation`: equation string with reactants and products separated by an equality operator.

# Returns

  - Tuple of (reactants_dict, products_dict, equal_sign_char).

# Examples

```julia
julia> reactants, products, sign = parse_equation("2H2 + O2 = 2H2O");

julia> reactants["H2"]
2

julia> products["H2O"]
2
```
"""
function parse_equation(equation::AbstractString)
    equal_sign = '='
    for c in equation
        if c in EQUAL_REACTION_SET
            equal_sign = c
            break
        end
    end

    sides = strip.(split(equation, EQUAL_REACTION))
    nsides = length(sides)
    left_side = nsides > 0 ? sides[1] : ""
    right_side = nsides > 1 ? sides[2] : ""

    function parse_side(side::AbstractString)
        terms = split(side, " +")
        result = OrderedDict{String,Real}()

        for term in terms
            t = strip(term)

            m = match(r"^(?<coeff>[-+]?\d+//\d+|[-+]?\d*\.?\d+)?\s*(?<formula>.+)$", t)

            if m !== nothing
                coeff_str = m[:coeff]
                formula = strip(m[:formula])

                coeff = if coeff_str === nothing || coeff_str == ""
                    1
                else
                    eval(Meta.parse(coeff_str))
                end

                if !(coeff isa Real)
                    error("Invalid coefficient: $coeff_str")
                end

                result[formula] = coeff
            else
                error("Unexpected term format: $term")
            end
        end

        return OrderedDict(k => stoich_coef_round(v) for (k, v) in result)
    end

    reactants = if left_side == "∅" || left_side == ""
        OrderedDict{String,Int}()
    else
        parse_side(left_side)
    end
    products = if right_side == "∅" || right_side == ""
        OrderedDict{String,Int}()
    else
        parse_side(right_side)
    end

    return reactants, products, equal_sign
end

"""
    colored_equation(equation::AbstractString) -> String

Generate a terminal-colored representation of a chemical equation.

# Arguments

  - `equation`: chemical equation string.

# Returns

  - String with ANSI color codes for terminal display.
"""
function colored_equation(equation::AbstractString)
    reactants, products, equal_sign = parse_equation(equation)
    left_side = if isempty(reactants)
        "∅"
    else
        join(
            [
                string(COL_STOICH_EXT(
                    if isone(v)
                        ""
                    elseif v < 0
                        "($(v))"
                    else
                        string(v)
                    end,
                )) * colored_formula(k) for (k, v) in reactants
            ],
            " + ",
        )
    end
    right_side = if isempty(products)
        "∅"
    else
        join(
            [
                string(COL_STOICH_EXT(
                    if isone(v)
                        ""
                    elseif v < 0
                        "($(v))"
                    else
                        string(v)
                    end,
                )) * colored_formula(k) for (k, v) in products
            ],
            " + ",
        )
    end
    return left_side * " " * string(COL_PAR(string(equal_sign))) * " " * right_side
end

"""
    format_equation(coeffs::AbstractDict; scaling=1, equal_sign='=') -> String

Format a stoichiometric coefficient dictionary into an equation string.

# Arguments

  - `coeffs`: dictionary mapping species to stoichiometric coefficients (negative = reactants, positive = products).
  - `scaling`: optional scaling factor for all coefficients (default 1).
  - `equal_sign`: equality operator to use (default '=').

# Returns

  - Formatted equation string with automatic electron balancing if needed.

# Examples

```julia
julia> coeffs = Dict("H2" => -2, "O2" => -1, "H2O" => 2);

julia> format_equation(coeffs)
"2H2 + O2 = 2H2O"
```    # Separate reactants and products
```
"""
function format_equation(coeffs::AbstractDict; scaling=1, equal_sign='=')
    # Separate reactants and products
    reactants = String[]
    products = String[]
    total_charge_left = 0
    total_charge_right = 0

    for (species, coeff) in coeffs
        if species !== "Zz"
            coeff = stoich_coef_round(coeff * scaling)

            # Format the coefficient
            abs_coeff = coeff < 0 ? -coeff : coeff
            coeff_str = if isapprox(abs_coeff, 1; atol=1e-6)
                ""
            elseif isinteger(abs_coeff)
                string(Int(abs_coeff))
            else
                string(abs_coeff)
            end

            # if coeff > 0
            #     push!(products, "$coeff_str$species")
            #     total_charge_right += coeff * extract_charge(species)
            # elseif coeff < 0
            #     push!(reactants, "$coeff_str$species")
            #     total_charge_left += coeff * extract_charge(species)
            # end

            if coeff > 0
                push!(products, "$coeff_str$species")
                total_charge_right += coeff * extract_charge(species)
            elseif coeff < 0
                push!(reactants, "$coeff_str$species")
                total_charge_left += coeff * extract_charge(species)
            elseif coeff == 0
            else
                push!(products, "($coeff_str)$species")
                total_charge_right += coeff * extract_charge(species)
            end
        end
    end

    # Build the initial equation
    left_side = join(reactants, " + ")
    right_side = join(products, " + ")

    # Compute the charge difference (corrected)
    charge_diff = total_charge_right + total_charge_left

    # Balance charges if necessary
    if !isapprox(charge_diff, 0; atol=1e-6)
        needed_e = stoich_coef_round(abs(charge_diff))
        e_term = needed_e == 1 ? "e⁻" : "$needed_e" * "e⁻"

        if charge_diff < 0
            # Add e- to the left (reactants)
            left_side = isempty(left_side) ? e_term : "$left_side + $e_term"
        else
            # Add e- to the right (products)
            right_side = isempty(right_side) ? e_term : "$right_side + $e_term"
        end
    end

    if length(left_side) == 0
        left_side = "∅"
    end
    if length(right_side) == 0
        right_side = "∅"
    end

    return "$left_side $(isnothing(equal_sign) ? '=' : equal_sign) $right_side"
end

"""
    remove_redundant_outer_parens_unicode(s::AbstractString) -> String

Remove unnecessary outer parentheses from a Unicode expression string.

# Arguments

  - `s`: expression string potentially wrapped in parentheses.

# Returns

  - String with redundant outer parentheses removed while preserving necessary ones.

Parentheses are kept if they contain root-level +/- operators or if removing them
would change the expression's meaning.
"""
function remove_redundant_outer_parens_unicode(s::AbstractString)
    while true
        if !startswith(s, "(") || !endswith(s, ")")
            break
        end

        # Use Unicode-safe indices
        first_paren = firstindex(s)
        last_paren = lastindex(s)

        # Find the index after first '(' and before last ')'
        inner_start = nextind(s, first_paren)
        inner_end = prevind(s, last_paren)

        # Check balanced parentheses over s
        count = 0
        balanced = true
        for c in s
            if c == '('
                count += 1
            elseif c == ')'
                count -= 1
                if count == 0 && c != last(s)  # parentheses close before end
                    balanced = false
                    break
                elseif count < 0
                    balanced = false
                    break
                end
            end
        end

        if !balanced || count != 0
            break
        end

        inner = s[inner_start:inner_end]

        # Check for '/' at root level (outside parentheses)
        function has_root_level_slash(str)
            lvl = 0
            for ch in str
                if ch == '('
                    lvl += 1
                elseif ch == ')'
                    lvl -= 1
                elseif ch == '/' && lvl == 0
                    return true
                end
            end
            return false
        end

        # Check for '+' or '-' at root level
        function has_root_level_plusminus(str)
            lvl = 0
            for ch in str
                if ch == '('
                    lvl += 1
                elseif ch == ')'
                    lvl -= 1
                elseif (ch == '+' || ch == '-') && lvl == 0
                    return true
                end
            end
            return false
        end

        if has_root_level_slash(inner)
            s = inner
        elseif has_root_level_plusminus(inner)
            break  # keep parentheses if + or - at root level
        else
            s = inner
        end
    end
    return s
end

"""
    add_parentheses_if_needed(s::String) -> String

Add parentheses to an expression string if it contains root-level +/- operators.

# Arguments

  - `s`: expression string.

# Returns

  - String wrapped in parentheses if needed, unchanged otherwise.

Parentheses are added only when the expression contains addition or subtraction
operators at the root level (outside any existing parentheses).
"""
function add_parentheses_if_needed(s::String)
    # Return early if s is already parenthesized fully and balanced
    if startswith(s, "(") && endswith(s, ")")
        # Check if outer parentheses fully enclose s
        count = 0
        for (i, c) in enumerate(s)
            if c == '('
                count += 1
            elseif c == ')'
                count -= 1
                if count == 0 && i != lastindex(s)
                    break
                end
            end
        end
        if count == 0
            return s  # already properly parenthesized
        end
    end

    # Helper to detect + or - at root level (outside any parentheses)
    function has_root_level_plusminus(str)
        lvl = 0
        for ch in str
            if ch == '('
                lvl += 1
            elseif ch == ')'
                lvl -= 1
            elseif (ch == '+' || ch == '-') && lvl == 0
                return true
            end
        end
        return false
    end

    # Add parentheses only if necessary
    if has_root_level_plusminus(s)
        return "(" * s * ")"
    else
        return s
    end
end
