"""
    ATOMIC_ORDER :: Vector{Symbol}

Canonical atomic ordering used across the package for serialization,
stringification and deterministic ordering of formula fields.

This vector lists element symbols in the preferred display/serialization
order and includes the special placeholder :Zz which represents a unit
positive charge in compositions.

# Examples

```julia
julia> print(ATOMIC_ORDER)
[:Ca, :Na, :K, :Mg, :Sr, :Ba, :Al, :Fe, :Ti, :Mn, :Cr, :Si, :C, :H, :N, :S, :O, :P, :B, :F, :Cl, :Br, :I, :Zz]
```
"""
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

"""
    CEMENT_TO_MENDELEEV :: Vector{Pair{Symbol,OrderedDict{Symbol,Int}}}

Mapping from cement shorthand symbols to their corresponding oxide or
elemental compositions. Each Pair maps a cement shorthand Symbol (key)
to an OrderedDict (value) describing composition in terms of element
symbols and integer stoichiometric coefficients.

This mapping is used to translate cement shorthand notation into full
elemental compositions for formula construction and serialization.

# Examples

```julia
julia> haskey(Dict(CEMENT_TO_MENDELEEV), :C)
true

julia> for (k,v) in CEMENT_TO_MENDELEEV println(k, " ≡ ", unicode(Species(v))) end
C ≡ CaO
M ≡ MgO
S ≡ SiO₂
A ≡ Al₂O₃
F ≡ Fe₂O₃
K ≡ K₂O
N ≡ Na₂O
P ≡ O₅P₂
T ≡ TiO₂
C̄ ≡ CO₂
S̄ ≡ SO₃
N̄ ≡ NO₃
H ≡ H₂O
```
"""
const CEMENT_TO_MENDELEEV = [
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

"""
    OXIDE_ORDER :: Vector{Symbol}

Derived ordered list of cement oxide shorthand symbols, extracted from
`CEMENT_TO_MENDELEEV` while preserving the original sequence. Useful for
deterministic iteration over oxide types.

# Examples

```julia
julia> print(OXIDE_ORDER)
[:C, :M, :S, :A, :F, :K, :N, :P, :T, :C̄, :S̄, :N̄, :H]
```
"""
const OXIDE_ORDER = collect(first.(CEMENT_TO_MENDELEEV))
