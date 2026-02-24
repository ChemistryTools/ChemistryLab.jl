"""
    struct ChemicalState{C, S, Q<:AbstractQuantity}

Immutable container holding the thermodynamic state of a `ChemicalSystem`.

Molar amounts are always stored internally in mol regardless of the input unit.
Each species can be provided independently as a molar amount (mol) or as a mass
(g, kg, etc.) — the constructor converts each entry individually using the
molar mass `M` stored in the corresponding species.

The struct itself is immutable — fields cannot be reassigned. However,
`n`, `T`, and `P` are stored as `Vector` to allow in-place mutation via
`set_moles!`, `set_temperature!`, and `set_pressure!`.

`system` is a shared reference: cloning via `Base.copy` does not duplicate
the underlying `ChemicalSystem`.

# Fields

  - `system`: reference to the underlying `ChemicalSystem`.
  - `n`: molar amounts (mol), one per species — mutable in place.
  - `T`: temperature (K) — 1-element Vector, mutable in place.
  - `P`: pressure (Pa) — 1-element Vector, mutable in place.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("Na+"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
       ]);

julia> state = ChemicalState(cs; T=298.15u"K", P=1u"bar");

julia> length(state.n)
2

julia> ustrip(state.T[])
298.15
```
"""
struct ChemicalState{C, S, Q<:AbstractQuantity}
    system::ChemicalSystem{C, S}    # shared reference — not duplicated on copy
    n::Vector{Q}                     # molar amounts [mol] — always stored in mol
    T::Vector{Q}                     # temperature [K]     — 1-element Vector for mutability
    P::Vector{Q}                     # pressure [Pa]       — 1-element Vector for mutability
end

# ── Internal helper ───────────────────────────────────────────────────────────

"""
    _entry_to_moles(v::AbstractQuantity, s::AbstractSpecies) -> AbstractQuantity

Convert a single value `v` to moles for species `s`.
If `v` has amount dimension (mol), it is returned as-is.
If `v` has mass dimension, it is divided by the molar mass `M` of `s`.
Otherwise an error is raised.
"""
function _entry_to_moles(v::AbstractQuantity, s::AbstractSpecies)
    if dimension(v) == dimension(u"mol")
        return v                            # already in mol — return as-is
    elseif dimension(v) == dimension(u"kg")
        return uconvert(us"mol", v / s[:M])  # m / M → mol, with unit conversion
    else
        error("Value for species $(symbol(s)) must have amount (mol) or mass dimension, got $(dimension(v))")
    end
end

# ── Constructors ──────────────────────────────────────────────────────────────

"""
    ChemicalState(system::ChemicalSystem; T, P, n) -> ChemicalState

Construct a `ChemicalState` from a `ChemicalSystem` with optional initial
temperature, pressure, and molar amounts (default: all zero).

# Arguments

  - `system`: the `ChemicalSystem` describing the species.
  - `T`: temperature in K (default: `298.15u"K"`).
  - `P`: pressure (default: `1u"bar"`).
  - `n`: molar amounts in mol (default: zeros).

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("Na+"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
       ]);

julia> state = ChemicalState(cs; T=298.15u"K", P=1u"bar");

julia> ustrip(state.T[])
298.15

julia> ustrip(state.P[]) ≈ 1e5
true

julia> all(iszero.(ustrip.(state.n)))
true
```
"""
function ChemicalState(
    system::ChemicalSystem;
    T::Q = 298.15u"K",
    P::Q = 1u"bar",
    n::AbstractVector = fill(0.0u"mol", length(system)),
) where {Q<:AbstractQuantity}
    @assert dimension(T) == dimension(u"K")  "T must have temperature dimension"
    @assert dimension(P) == dimension(u"Pa") "P must have pressure dimension"
    @assert length(n) == length(system) "n must have one entry per species (got $(length(n)), expected $(length(system)))"

    # Convert each entry individually — species i may be in mol while species j is in g
    n_mol = Q[uconvert(us"mol", _entry_to_moles(nᵢ, s))
              for (nᵢ, s) in zip(n, system.species)]

    return ChemicalState(system, n_mol, Q[T], Q[P])
end

"""
    ChemicalState(system::ChemicalSystem, values::AbstractVector; T, P) -> ChemicalState

Construct a `ChemicalState` with explicit initial amounts or masses.

Each entry of `values` is converted to moles independently: it can be a molar
amount (mol) or a mass (g, kg, etc.). Different entries may use different units.

# Arguments

  - `system`: the `ChemicalSystem`.
  - `values`: one value per species — mol or mass, mixed units allowed.
  - `T`: temperature (default: `298.15u"K"`).
  - `P`: pressure (default: `1u"bar"`).

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("NaCl"; aggregate_state=AS_CRYSTAL),
       ]);

julia> state = ChemicalState(cs, [55.5u"mol", 5.844u"g"]);

julia> ustrip(moles(state, "H2O"))
55.5

julia> isapprox(ustrip(moles(state, "NaCl")), 0.1; rtol=1.e-4)
true
```
"""
function ChemicalState(
    system::ChemicalSystem,
    values::AbstractVector;
    T::Q = 298.15u"K",
    P::Q = 1u"bar",
) where {Q<:AbstractQuantity}
    return ChemicalState(system; T=T, P=P, n=values)
end

# ── Temperature and pressure accessors ───────────────────────────────────────

"""
    temperature(state::ChemicalState) -> AbstractQuantity

Return the current temperature.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> state = ChemicalState(cs; T=298.15u"K", P=1u"bar");

julia> ustrip(temperature(state))
298.15
```
"""
temperature(state::ChemicalState) = state.T[]  # dereference the 1-element Vector

"""
    pressure(state::ChemicalState) -> AbstractQuantity

Return the current pressure.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> state = ChemicalState(cs; T=298.15u"K", P=1u"bar");

julia> isapprox(ustrip(pressure(state)), 1e5; rtol=1.e-4)
true
```
"""
pressure(state::ChemicalState) = state.P[]  # dereference the 1-element Vector

"""
    set_temperature!(state::ChemicalState, T::AbstractQuantity) -> ChemicalState

Set the temperature in place.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> state = ChemicalState(cs; T=298.15u"K", P=1u"bar");

julia> set_temperature!(state, 350.0u"K");

julia> ustrip(temperature(state))
350.0
```
"""
function set_temperature!(state::ChemicalState, T::AbstractQuantity)
    @assert dimension(T) == dimension(u"K") "T must have temperature dimension"
    state.T[] = T   # mutate the 1-element Vector in place
    return state
end

"""
    set_pressure!(state::ChemicalState, P::AbstractQuantity) -> ChemicalState

Set the pressure in place.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> state = ChemicalState(cs; T=298.15u"K", P=1u"bar");

julia> set_pressure!(state, 2u"bar");

julia> isapprox(ustrip(pressure(state)), 2e5; rtol=1.e-4)
true
```
"""
function set_pressure!(state::ChemicalState, P::AbstractQuantity)
    @assert dimension(P) == dimension(u"Pa") "P must have pressure dimension"
    state.P[] = P   # mutate the 1-element Vector in place
    return state
end

# ── Molar amount accessors ────────────────────────────────────────────────────

"""
    moles(state::ChemicalState, s::AbstractSpecies) -> AbstractQuantity

Return the molar amount of species `s` in mol.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> state = ChemicalState(cs, [55.5u"mol"]);

julia> ustrip(moles(state, cs[1]))
55.5
```
"""
function moles(state::ChemicalState, s::AbstractSpecies)
    i = findfirst(x -> x == s, state.system.species)   # lookup index by equality
    isnothing(i) && error("Species $(symbol(s)) not found in ChemicalSystem")
    return state.n[i]
end

"""
    moles(state::ChemicalState, sym::AbstractString) -> AbstractQuantity

Return the molar amount of the species identified by symbol `sym` in mol.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> state = ChemicalState(cs, [55.5u"mol"]);

julia> ustrip(moles(state, "H2O"))
55.5
```
"""
function moles(state::ChemicalState, sym::AbstractString)
    return moles(state, state.system[sym])  # resolve symbol via dict then delegate
end

"""
    mass(state::ChemicalState, s::AbstractSpecies) -> AbstractQuantity

Return the mass of species `s`, computed as `n × M`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> state = ChemicalState(cs, [55.5u"mol"]);

julia> ustrip(uconvert(us"g", mass(state, cs[1]))) ≈ 55.5 * 18.015
true
```
"""
function mass(state::ChemicalState, s::AbstractSpecies)
    i = findfirst(x -> x == s, state.system.species)
    isnothing(i) && error("Species $(symbol(s)) not found in ChemicalSystem")
    return state.n[i] * s[:M]   # n [mol] × M [g/mol] → mass
end

"""
    mass(state::ChemicalState, sym::AbstractString) -> AbstractQuantity

Return the mass of the species identified by symbol `sym`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> state = ChemicalState(cs, [55.5u"mol"]);

julia> ustrip(uconvert(us"g", mass(state, "H2O"))) ≈ 55.5 * 18.015
true
```
"""
function mass(state::ChemicalState, sym::AbstractString)
    return mass(state, state.system[sym])   # resolve symbol then delegate
end

"""
    set_moles!(state::ChemicalState, s::AbstractSpecies, n::AbstractQuantity) -> ChemicalState

Set the molar amount of species `s` in place.
If `n` has mass dimension, it is automatically converted to moles using `M`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> state = ChemicalState(cs, [55.5u"mol"]);

julia> set_moles!(state, cs[1], 10.0u"mol");

julia> ustrip(moles(state, "H2O"))
10.0
```
"""
function set_moles!(state::ChemicalState, s::AbstractSpecies, n::AbstractQuantity)
    i = findfirst(x -> x == s, state.system.species)
    isnothing(i) && error("Species $(symbol(s)) not found in ChemicalSystem")
    # Convert to mol in place — accepts mol or any mass unit
    state.n[i] = uconvert(us"mol", _entry_to_moles(n, s))
    return state
end

"""
    set_moles!(state::ChemicalState, sym::AbstractString, n::AbstractQuantity) -> ChemicalState

Set the molar amount of the species identified by symbol `sym` in place.
If `n` has mass dimension, it is automatically converted to moles using `M`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> state = ChemicalState(cs, [55.5u"mol"]);

julia> set_moles!(state, "H2O", 10.0u"mol");

julia> ustrip(moles(state, "H2O"))
10.0
```
"""
function set_moles!(state::ChemicalState, sym::AbstractString, n::AbstractQuantity)
    return set_moles!(state, state.system[sym], n)  # resolve symbol then delegate
end

# ── Clone ─────────────────────────────────────────────────────────────────────

"""
    Base.copy(state::ChemicalState) -> ChemicalState

Create a clone of a `ChemicalState` that shares the same `ChemicalSystem`
reference but owns independent copies of `n`, `T`, and `P`.

Modifying the clone does not affect the original, and vice versa.
The underlying `ChemicalSystem` is not duplicated.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("Na+"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
       ]);

julia> state = ChemicalState(cs, [55.5u"mol", 0.1u"mol"]);

julia> clone = copy(state);

julia> set_moles!(clone, "Na+", 0.5u"mol");

julia> ustrip(moles(state, "Na+"))
0.1

julia> ustrip(moles(clone, "Na+"))
0.5

julia> clone.system === state.system
true
```
"""
function Base.copy(state::ChemicalState)
    return ChemicalState(
        state.system,       # shared reference — ChemicalSystem is not duplicated
        copy(state.n),      # independent copy of molar amounts
        copy(state.T),      # independent copy of temperature
        copy(state.P),      # independent copy of pressure
    )
end

# ── Display ───────────────────────────────────────────────────────────────────

"""
    Base.show(io::IO, ::MIME"text/plain", state::ChemicalState)

Detailed multi-line REPL display for `ChemicalState`.
Shows both molar amounts and masses for each species.
"""
function Base.show(io::IO, ::MIME"text/plain", state::ChemicalState)
    pad = maximum(length.(unicode.(state.system.species)); init=8)
    println(io, typeof(state))
    println(io, lpad("T", pad), " : ", state.T[])
    println(io, lpad("P", pad), " : ", uconvert(us"bar", state.P[]))
    println(io, lpad("species", pad), " : ", rpad("n [mol]", 25), "m [g]")
    println(io, repeat("─", pad + 50))
    for (s, nᵢ) in zip(state.system.species, state.n)
        mᵢ = uconvert(us"g", nᵢ * s[:M])    # compute mass from moles and molar mass
        println(io, lpad(unicode(s), pad), " : ", rpad(string(nᵢ), 25), mᵢ)
    end
end

"""
    Base.show(io::IO, state::ChemicalState)

Compact single-line representation of a `ChemicalState`.
"""
function Base.show(io::IO, state::ChemicalState)
    return print(io, "ChemicalState(T=", state.T[], ", P=", state.P[],
                 ", ", length(state.n), " species)")
end
