"""
    struct ChemicalState{C, S, Q<:AbstractQuantity}

Immutable container holding the thermodynamic state of a `ChemicalSystem`.

Molar amounts are always stored internally in mol regardless of the input unit.
Each species can be provided independently as a molar amount (mol) or as a mass
(g, kg, etc.) — the constructor converts each entry individually using the
molar mass `M` stored in the corresponding species.

The struct itself is immutable — fields cannot be reassigned. However,
`n`, `T`, and `P` are stored as `Vector` to allow in-place mutation via
`set_quantity!`, `set_temperature!`, and `set_pressure!`.

`system` is a shared reference: cloning via `Base.copy` does not duplicate
the underlying `ChemicalSystem`.

# Fields

  - `system`: reference to the underlying `ChemicalSystem`.
  - `n`: molar amounts (mol), one per species — mutable in place.
  - `T`: temperature (K) — 1-element Vector, mutable in place.
  - `P`: pressure (Pa) — 1-element Vector, mutable in place.
  - `n_phases`: moles per phase `(liquid, solid, gas, total)`.
  - `m_phases`: mass per phase `(liquid, solid, gas, total)`.
  - `V_phases`: volume per phase `(liquid, solid, gas, total)`.
  - `pH`: pH of the liquid phase, or `nothing` if H⁺ is absent.
  - `pOH`: pOH of the liquid phase, or `nothing` if OH⁻ is absent.
  - `porosity`: `(V_liquid + V_gas) / V_total`, or `nothing` if volumes unavailable.
  - `saturation`: `V_liquid / (V_liquid + V_gas)`, or `nothing` if pore volume is zero.

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
    system::ChemicalSystem{C, S}        # shared reference — not duplicated on copy
    n::Vector{Q}                         # molar amounts [mol] — always stored in mol
    T::Vector{Q}                         # temperature [K]     — 1-element Vector for mutability
    P::Vector{Q}                         # pressure [Pa]       — 1-element Vector for mutability
    n_phases::Vector{NamedTuple}         # moles per phase — 1-element Vector for mutability
    m_phases::Vector{NamedTuple}         # mass per phase  — 1-element Vector for mutability
    V_phases::Vector{NamedTuple}         # volume per phase — 1-element Vector for mutability
    pH::Vector{Union{Number, Nothing}}          # dimensionless — Float64 or Dual or nothing
    pOH::Vector{Union{Number, Nothing}}         # dimensionless — Float64 or Dual or nothing
    porosity::Vector{Union{Number, Nothing}}    # dimensionless ratio — V_pore / V_total
    saturation::Vector{Union{Number, Nothing}}  # dimensionless ratio — V_liq / V_pore
end

# ── Internal helpers ──────────────────────────────────────────────────────────

"""
    _entry_to_moles(v::AbstractQuantity, s::AbstractSpecies) -> AbstractQuantity

Convert a single value `v` to moles for species `s`.
If `v` has amount dimension (mol), it is returned as-is.
If `v` has mass dimension, it is divided by the molar mass `M` of `s`.
Otherwise an error is raised.
"""
function _entry_to_moles(v::AbstractQuantity, s::AbstractSpecies)
    if dimension(v) == dimension(u"mol")
        return v                                # already in mol — return as-is
    elseif dimension(v) == dimension(u"kg")
        return uconvert(us"mol", v / s[:M])     # m / M → mol, with unit conversion
    else
        error("Value for species $(symbol(s)) must have amount (mol) or mass dimension, got $(dimension(v))")
    end
end

"""
    _has_molar_volume(s::AbstractSpecies) -> Bool

Return `true` if species `s` has a standard molar volume `V⁰` available.
"""
_has_molar_volume(s::AbstractSpecies) = haskey(s, :V⁰)

"""
    _molar_volume(s::AbstractSpecies) -> ThermoFunction

Return the standard molar volume ThermoFunction of species `s`.
Must be called as `_molar_volume(s)(T=T, P=P; unit=true)` to get a quantity.
"""
_molar_volume(s::AbstractSpecies) = s[:V⁰]

"""
    _compute_n_phases(system, n) -> NamedTuple

Compute moles per phase from species vector `n`.
"""
function _compute_n_phases(system::ChemicalSystem, n::AbstractVector)
    function _phase(idx)
        return sum((n[i] for i in idx); init=0.0u"mol")
    end
    n_liquid = _phase(system.idx_aqueous)
    n_solid  = _phase(system.idx_crystal)
    n_gas    = _phase(system.idx_gas)
    return (liquid=n_liquid, solid=n_solid, gas=n_gas, total=n_liquid + n_solid + n_gas)
end

"""
    _compute_m_phases(system, n) -> NamedTuple

Compute mass per phase from species vector `n` and molar masses.
"""
function _compute_m_phases(system::ChemicalSystem, n::AbstractVector)
    function _phase(idx)
        return sum((n[i] * system.species[i][:M] for i in idx); init=0.0u"kg")
    end
    m_liquid = _phase(system.idx_aqueous)
    m_solid  = _phase(system.idx_crystal)
    m_gas    = _phase(system.idx_gas)
    return (liquid=m_liquid, solid=m_solid, gas=m_gas, total=m_liquid + m_solid + m_gas)
end

"""
    _compute_V_phases(system, n, T, P) -> NamedTuple

Compute volume per phase from `n`, `T`, `P` and standard molar volumes `V⁰`.
Gas phase falls back to ideal gas law if `V⁰` is not available for all gas species.
"""
function _compute_V_phases(system::ChemicalSystem, n::AbstractVector, T, P)
    # Sum n × V⁰(T,P) over indices where V⁰ is available
    _phase(idx) = sum(
        n[i] * _molar_volume(system.species[i])(T=T, P=P; unit=true)
        for i in idx if _has_molar_volume(system.species[i]);
        init=0.0u"m^3",
    )

    V_liquid = _phase(system.idx_aqueous)
    V_solid  = _phase(system.idx_crystal)

    if isempty(system.idx_gas)
        V_gas = 0.0u"m^3"
    elseif all(_has_molar_volume(system.species[i]) for i in system.idx_gas)
        V_gas = _phase(system.idx_gas)          # use database values if all available
    else
        R     = Constants.R                     # ideal gas constant
        n_gas = sum(n[i] for i in system.idx_gas; init=0.0u"mol")
        V_gas = uconvert(u"m^3", n_gas * R * T / P)    # ideal gas fallback
    end

    V_total = V_liquid + V_solid + V_gas
    return (liquid=V_liquid, solid=V_solid, gas=V_gas, total=V_total)
end

"""
    _compute_porosity(V_phases) -> Union{Float64, Nothing}

Compute porosity = (V_liquid + V_gas) / V_total.
Returns `nothing` if total volume is zero.
"""
function _compute_porosity(V_phases)
    V_tot  = V_phases.total
    iszero(ustrip(V_tot)) && return nothing
    return (V_phases.liquid + V_phases.gas) / V_tot   # dimensionless
end

"""
    _compute_saturation(V_phases) -> Union{Float64, Nothing}

Compute saturation = V_liquid / (V_liquid + V_gas).
Returns `nothing` if pore volume is zero.
"""
function _compute_saturation(V_phases)
    V_pore = V_phases.liquid + V_phases.gas
    iszero(ustrip(V_pore)) && return nothing
    return V_phases.liquid / V_pore                   # dimensionless
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

    # Compute all derived quantities once at construction time
    n_ph = _compute_n_phases(system, n_mol)
    m_ph = _compute_m_phases(system, n_mol)
    V_ph = _compute_V_phases(system, n_mol, T, P)
    _pH = _compute_pH( system, n_mol, T, P, V_ph.liquid)
    _pOH = _compute_pOH(system, n_mol, T, P, V_ph.liquid)
    _porosity = _compute_porosity(V_ph)
    _saturation = _compute_saturation(V_ph)

    return ChemicalState(
        system,
        n_mol,
        Q[T],
        Q[P],
        NamedTuple[n_ph],       # 1-element Vector for in-place update
        NamedTuple[m_ph],
        NamedTuple[V_ph],
        Union{Number, Nothing}[_pH],
        Union{Number, Nothing}[_pOH],
        Union{Number, Nothing}[_porosity],
        Union{Number, Nothing}[_saturation],
    )
end

"""
    ChemicalState(system::ChemicalSystem, values::AbstractVector; T, P) -> ChemicalState

Construct a `ChemicalState` with explicit initial amounts or masses.
Each entry is converted to moles independently — mixed units allowed.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("NaCl"; aggregate_state=AS_CRYSTAL),
       ]);

julia> state = ChemicalState(cs, [55.5u"mol", 5.844u"g"]);

julia> ustrip(moles(state, "H2O"))
55.5

julia> isapprox(ustrip(moles(state, "NaCl")), 0.1; rtol=1e-4)
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

# ── Internal update ───────────────────────────────────────────────────────────

"""
    _update_derived!(state::ChemicalState)

Recompute and update in place all derived quantities after any mutation
of `n`, `T`, or `P`. Called automatically by `set_quantity!`,
`set_temperature!`, and `set_pressure!`.
"""
function _update_derived!(state::ChemicalState)
    T    = temperature(state)
    P    = pressure(state)
    V_ph = _compute_V_phases(state.system, state.n, T, P)
    state.n_phases[]   = _compute_n_phases(state.system, state.n)
    state.m_phases[]   = _compute_m_phases(state.system, state.n)
    state.V_phases[]   = V_ph
    state.pH[]         = _compute_pH( state.system, state.n, T, P, V_ph.liquid)
    state.pOH[]        = _compute_pOH(state.system, state.n, T, P, V_ph.liquid)
    state.porosity[]   = _compute_porosity(V_ph)
    state.saturation[] = _compute_saturation(V_ph)
    return state
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
temperature(state::ChemicalState) = state.T[]

"""
    pressure(state::ChemicalState) -> AbstractQuantity

Return the current pressure.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> state = ChemicalState(cs; T=298.15u"K", P=1u"bar");

julia> isapprox(ustrip(pressure(state)), 1e5; rtol=1e-4)
true
```
"""
pressure(state::ChemicalState) = state.P[]

"""
    set_temperature!(state::ChemicalState, T::AbstractQuantity) -> ChemicalState

Set the temperature in place and update all derived quantities.

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
    state.T[] = T
    _update_derived!(state)     # volumes depend on T — recompute everything
    return state
end

"""
    set_pressure!(state::ChemicalState, P::AbstractQuantity) -> ChemicalState

Set the pressure in place and update all derived quantities.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> state = ChemicalState(cs; T=298.15u"K", P=1u"bar");

julia> set_pressure!(state, 2u"bar");

julia> isapprox(ustrip(pressure(state)), 2e5; rtol=1e-4)
true
```
"""
function set_pressure!(state::ChemicalState, P::AbstractQuantity)
    @assert dimension(P) == dimension(u"Pa") "P must have pressure dimension"
    state.P[] = P
    _update_derived!(state)     # volumes depend on P — recompute everything
    return state
end

# ── Molar amount accessors ────────────────────────────────────────────────────

"""
    moles(state::ChemicalState) -> NamedTuple

Return moles per phase `(liquid, solid, gas, total)`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O";  aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("NaCl"; aggregate_state=AS_CRYSTAL),
       ]);

julia> state = ChemicalState(cs, [55.5u"mol", 0.05u"mol"]);

julia> ustrip(moles(state).liquid)
55.5
```
"""
moles(state::ChemicalState) = state.n_phases[]

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
    i = findfirst(x -> x == s, state.system.species)
    isnothing(i) && error("Species $(symbol(s)) not found in ChemicalSystem")
    return state.n[i]
end

"""
    moles(state::ChemicalState, sym::AbstractString) -> AbstractQuantity

Return the molar amount of the species identified by symbol `sym`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> state = ChemicalState(cs, [55.5u"mol"]);

julia> ustrip(moles(state, "H2O"))
55.5
```
"""
moles(state::ChemicalState, sym::AbstractString) = moles(state, state.system[sym])

# ── Mass accessors ────────────────────────────────────────────────────────────

"""
    mass(state::ChemicalState) -> NamedTuple

Return mass per phase `(liquid, solid, gas, total)`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O";  aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("NaCl"; aggregate_state=AS_CRYSTAL),
       ]);

julia> state = ChemicalState(cs, [55.5u"mol", 0.05u"mol"]);

julia> mass(state).total isa AbstractQuantity
true
```
"""
mass(state::ChemicalState) = state.m_phases[]

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
    return state.n[i] * s[:M]
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
mass(state::ChemicalState, sym::AbstractString) = mass(state, state.system[sym])

# ── Volume accessors ──────────────────────────────────────────────────────────

"""
    volume(state::ChemicalState) -> NamedTuple

Return volume per phase `(liquid, solid, gas, total)`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O";  aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("NaCl"; aggregate_state=AS_CRYSTAL),
       ]);

julia> state = ChemicalState(cs, [55.5u"mol", 0.05u"mol"]);

julia> volume(state).total isa AbstractQuantity
true
```
"""
volume(state::ChemicalState) = state.V_phases[]

"""
    volume(state::ChemicalState, s::AbstractSpecies) -> Union{AbstractQuantity, Nothing}

Return the volume contribution of species `s` as `n × V⁰(T,P)`.
Returns `nothing` if `V⁰` is not available for `s`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> state = ChemicalState(cs, [55.5u"mol"]);

julia> volume(state, cs[1]) isa Union{AbstractQuantity, Nothing}
true
```
"""
function volume(state::ChemicalState, s::AbstractSpecies)
    i = findfirst(x -> x == s, state.system.species)
    isnothing(i) && error("Species $(symbol(s)) not found in ChemicalSystem")
    _has_molar_volume(s) || return nothing
    return state.n[i] * _molar_volume(s)(T=temperature(state), P=pressure(state); unit=true)
end

"""
    volume(state::ChemicalState, sym::AbstractString) -> Union{AbstractQuantity, Nothing}

Return the volume contribution of the species identified by symbol `sym`.
Returns `nothing` if `V⁰` is not available.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> state = ChemicalState(cs, [55.5u"mol"]);

julia> volume(state, "H2O") isa Union{AbstractQuantity, Nothing}
true
```
"""
volume(state::ChemicalState, sym::AbstractString) = volume(state, state.system[sym])

# ── pH, pOH, porosity, saturation accessors ───────────────────────────────────

"""
    pH(state::ChemicalState) -> Union{Float64, Nothing}

Return the pH of the liquid phase, or `nothing` if H⁺ is absent.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("H+";  aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
       ]);

julia> state = ChemicalState(cs, [55.5u"mol", 1e-7u"mol"]);

julia> pH(state) isa Union{Float64, Nothing}
true
```
"""
pH(state::ChemicalState) = state.pH[]

"""
    pOH(state::ChemicalState) -> Union{Float64, Nothing}

Return the pOH of the liquid phase, or `nothing` if OH⁻ is absent.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("OH-"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
       ]);

julia> state = ChemicalState(cs, [55.5u"mol", 1e-7u"mol"]);

julia> pOH(state) isa Union{Float64, Nothing}
true
```
"""
pOH(state::ChemicalState) = state.pOH[]

"""
    porosity(state::ChemicalState) -> Union{Float64, Nothing}

Return the porosity `(V_liquid + V_gas) / V_total`,
or `nothing` if total volume is zero.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O";  aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("NaCl"; aggregate_state=AS_CRYSTAL),
       ]);

julia> state = ChemicalState(cs, [55.5u"mol", 0.05u"mol"]);

julia> porosity(state) isa Union{Float64, Nothing}
true
```
"""
porosity(state::ChemicalState) = state.porosity[]

"""
    saturation(state::ChemicalState) -> Union{Float64, Nothing}

Return the saturation `V_liquid / (V_liquid + V_gas)`,
or `nothing` if pore volume is zero.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O";  aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("NaCl"; aggregate_state=AS_CRYSTAL),
       ]);

julia> state = ChemicalState(cs, [55.5u"mol", 0.05u"mol"]);

julia> saturation(state) isa Union{Float64, Nothing}
true
```
"""
saturation(state::ChemicalState) = state.saturation[]

"""
    _compute_pH(system, n, V_liquid) -> Union{Float64, Nothing}

Compute pH using the most reliable species between H⁺ and OH⁻.

- If c_H⁺ ≥ c_OH⁻: pH = -log10(c_H⁺)
- If c_OH⁻ > c_H⁺: pOH = -log10(c_OH⁻), then pH = pKw(T) - pOH

pKw(T) is retrieved from the reaction H2O = H+ + OH- in the system's
reaction dictionary, evaluated at the current T and P.
Returns `nothing` if neither H⁺ nor OH⁻ is present, or if liquid volume is zero.
"""
function _compute_pH(system::ChemicalSystem, n::AbstractVector, T, P, V_liquid)
    iszero(ustrip(V_liquid)) && return nothing

    i_H  = findfirst(s -> symbol(s) == "H+",  system.species)
    i_OH = findfirst(s -> symbol(s) == "OH-", system.species)

    # Need at least one of H+ or OH-
    isnothing(i_H) && isnothing(i_OH) && return nothing

    c_H  = isnothing(i_H)  ? 0.0 : ustrip(us"mol/L", n[i_H]  / V_liquid)
    c_OH = isnothing(i_OH) ? 0.0 : ustrip(us"mol/L", n[i_OH] / V_liquid)

    if c_H >= c_OH
        # Acidic or neutral — pH directly from H+
        c_H <= 0.0 && return nothing
        return -log10(c_H)
    else
        # Basic — use OH- and pKw
        c_OH <= 0.0 && return nothing
        pOH = -log10(c_OH)
        pKw = _compute_pKw(system, T, P)
        isnothing(pKw) && return pOH   # fallback: cannot reconstruct pH without pKw
        return pKw - pOH
    end
end

"""
    _compute_pOH(system, n, T, P, V_liquid) -> Union{Float64, Nothing}

Compute pOH symmetrically to `_compute_pH`:

- If c_OH⁻ ≥ c_H⁺: pOH = -log10(c_OH⁻)
- If c_H⁺ > c_OH⁻: pH = -log10(c_H⁺), then pOH = pKw(T) - pH

Returns `nothing` if neither species is present or volume is zero.
"""
function _compute_pOH(system::ChemicalSystem, n::AbstractVector, T, P, V_liquid)
    iszero(ustrip(V_liquid)) && return nothing

    i_H  = findfirst(s -> symbol(s) == "H+",  system.species)
    i_OH = findfirst(s -> symbol(s) == "OH-", system.species)

    isnothing(i_H) && isnothing(i_OH) && return nothing

    c_H  = isnothing(i_H)  ? 0.0 : ustrip(us"mol/L", n[i_H]  / V_liquid)
    c_OH = isnothing(i_OH) ? 0.0 : ustrip(us"mol/L", n[i_OH] / V_liquid)

    if c_OH >= c_H
        c_OH <= 0.0 && return nothing
        return -log10(c_OH)
    else
        c_H <= 0.0 && return nothing
        pH  = -log10(c_H)
        pKw = _compute_pKw(system, T, P)
        isnothing(pKw) && return pH    # fallback
        return pKw - pH
    end
end

"""
    _compute_pKw(system::ChemicalSystem, T, P) -> Union{Float64, Nothing}

Compute pKw = -logK⁰(T, P) for the water dissociation reaction
H2O@ = H+ + OH- reconstructed on the fly from the species present in `system`.

Returns `nothing` if any of H2O@, H+, or OH- is absent from the system.
"""
function _compute_pKw(system::ChemicalSystem, T, P)
    # All three species must be present
    i_H2O = findfirst(s -> symbol(s) == "H2O@", system.species)
    i_H   = findfirst(s -> symbol(s) == "H+",   system.species)
    i_OH  = findfirst(s -> symbol(s) == "OH-",  system.species)

    (isnothing(i_H2O) || isnothing(i_H) || isnothing(i_OH)) && return nothing

    # Reconstruct the reaction H2O@ → H+ + OH- from the species objects
    r = system.species[i_H2O] → system.species[i_H] + system.species[i_OH]

    return -r.logK⁰(T=T, P=P)
end

# ── Mutation ──────────────────────────────────────────────────────────────────

"""
    set_quantity!(state::ChemicalState, s::AbstractSpecies, n::AbstractQuantity) -> ChemicalState

Set the molar amount of species `s` in place and update all derived quantities.
If `n` has mass dimension, it is automatically converted to moles using `M`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> state = ChemicalState(cs, [55.5u"mol"]);

julia> set_quantity!(state, cs[1], 10.0u"mol");

julia> ustrip(moles(state, "H2O"))
10.0
```
"""
function set_quantity!(state::ChemicalState, s::AbstractSpecies, n::AbstractQuantity)
    i = findfirst(x -> x == s, state.system.species)
    isnothing(i) && error("Species $(symbol(s)) not found in ChemicalSystem")
    state.n[i] = uconvert(us"mol", _entry_to_moles(n, s))
    _update_derived!(state)     # recompute all derived quantities
    return state
end

"""
    set_quantity!(state::ChemicalState, sym::AbstractString, n::AbstractQuantity) -> ChemicalState

Set the molar amount of the species identified by symbol `sym` in place
and update all derived quantities.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> state = ChemicalState(cs, [55.5u"mol"]);

julia> set_quantity!(state, "H2O", 10.0u"mol");

julia> ustrip(moles(state, "H2O"))
10.0
```
"""
set_quantity!(state::ChemicalState, sym::AbstractString, n::AbstractQuantity) =
    set_quantity!(state, state.system[sym], n)

# ── Clone ─────────────────────────────────────────────────────────────────────

"""
    Base.copy(state::ChemicalState) -> ChemicalState

Create a clone of a `ChemicalState` that shares the same `ChemicalSystem`
reference but owns independent copies of all mutable fields.

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

julia> set_quantity!(clone, "Na+", 0.5u"mol");

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
        state.system,                               # shared reference — not duplicated
        copy(state.n),                              # independent copy of molar amounts
        copy(state.T),                              # independent copy of temperature
        copy(state.P),                              # independent copy of pressure
        NamedTuple[state.n_phases[]],               # independent copy of phase moles
        NamedTuple[state.m_phases[]],               # independent copy of phase masses
        NamedTuple[state.V_phases[]],               # independent copy of phase volumes
        Union{Number, Nothing}[state.pH[]],         # independent copy of pH
        Union{Number, Nothing}[state.pOH[]],        # independent copy of pOH
        Union{Number, Nothing}[state.porosity[]],   # independent copy of porosity
        Union{Number, Nothing}[state.saturation[]], # independent copy of saturation
    )
end

# ── Display ───────────────────────────────────────────────────────────────────

"""
    Base.show(io::IO, ::MIME"text/plain", state::ChemicalState)

Detailed multi-line display for `ChemicalState`.
Shows molar amounts, masses, and volumes for each species grouped by phase,
with phase totals and scalar diagnostics (pH, pOH, porosity, saturation).
"""
function Base.show(io::IO, ::MIME"text/plain", state::ChemicalState)
    cs  = state.system
    T   = temperature(state)
    P   = pressure(state)

    # pad must accommodate species symbols AND phase labels
    phase_labels = ["# liquid #", "tot. liquid", "# solid #", "tot. solid",
                    "# gas #",    "tot. gas",    "# TOTAL #"]
    pad = maximum(length.([symbol.(cs.species); phase_labels]); init=12)

    show_volume = any(_has_molar_volume(s) for s in cs.species)
    show_conc   = show_volume && !isempty(cs.idx_aqueous)   # needs V_liquid
    show_ppart  = !isempty(cs.idx_gas)                      # needs only n_gas, P

    col_n = 20
    col_m = 20
    col_v = 20
    col_c = 20   # c [mol/L] — liquid only
    col_p = 20   # p [bar]   — gas only

    total_width = 1 + pad + 1 + col_n + 1 + col_m +
                  (show_volume ? 1 + col_v : 0) +
                  (show_conc   ? 1 + col_c : 0) +
                  (show_ppart  ? 1 + col_p : 0) + 1

    hl  = "─"; hhl = "═"; vl  = "│"
    tl  = "┌"; tr  = "┐"; bl  = "└"; br  = "┘"
    ml  = "├"; mr  = "┤"; mml = "╞"; mmr = "╡"
    ht  = "┄"

    _hline()   = ml  * repeat(hl,  total_width - 2) * mr
    _htline()  = ml  * repeat(ht,  total_width - 2) * mr
    _hhline()  = mml * repeat(hhl, total_width - 2) * mmr
    _topline() = tl  * repeat(hl,  total_width - 2) * tr
    _botline() = bl  * repeat(hl,  total_width - 2) * br

    _rpad(s, n) = s * repeat(" ", max(0, n - length(s)))
    _lpad(s, n) = repeat(" ", max(0, n - length(s))) * s

    _fmt(x)  = string(round(Float64(x); sigdigits=6))
    _fmt4(x) = string(round(Float64(x); digits=4))
    _fmt6(x) = string(round(Float64(x); digits=6))

    function _row(col1, cell_n, cell_m, cell_v, cell_c, cell_p)
        r = vl * _lpad(col1, pad) *
            vl * _lpad(cell_n, col_n) *
            vl * _lpad(cell_m, col_m)
        show_volume && (r *= vl * _lpad(cell_v, col_v))
        show_conc   && (r *= vl * _lpad(cell_c, col_c))
        show_ppart  && (r *= vl * _lpad(cell_p, col_p))
        r *= vl
        return r
    end

    inner = total_width - 2
    function _full_row(label, val)
        content = _lpad(label, pad) * " : " * val
        return vl * _rpad(content, inner) * vl
    end

    # ── Header ────────────────────────────────────────────────────────────────
    println(io, typeof(state))
    println(io, _topline())
    println(io, _full_row("T", _fmt4(ustrip(us"K",   T)) * " K"))
    println(io, _full_row("P", _fmt4(ustrip(us"bar", P)) * " bar"))

    # ── Helpers ───────────────────────────────────────────────────────────────

    V_liq = state.V_phases[].liquid
    function _concentration(nᵢ)
        !show_conc                && return ""
        iszero(ustrip(V_liq))    && return "N/A"
        return _fmt(ustrip(us"mol/L", nᵢ / V_liq))
    end

    n_gas_total = state.n_phases[].gas
    function _partial_pressure(nᵢ)
        !show_ppart                     && return ""
        iszero(ustrip(n_gas_total))     && return "N/A"
        xᵢ = ustrip(us"mol", nᵢ) / ustrip(us"mol", n_gas_total)
        return _fmt(xᵢ * ustrip(us"bar", P))
    end

    function _species_row(s, nᵢ, phase_key)
        n_val = _fmt(ustrip(us"mol", nᵢ))
        m_val = _fmt(ustrip(us"g",   nᵢ * s[:M]))
        V_val = if show_volume && _has_molar_volume(s)
            _fmt(ustrip(us"cm^3", nᵢ * _molar_volume(s)(T=T, P=P; unit=true)))
        elseif show_volume
            "N/A"
        else "" end
        c_val = phase_key == :liquid ? _concentration(nᵢ) : ""
        p_val = phase_key == :gas    ? _partial_pressure(nᵢ) : ""
        println(io, _row(symbol(s), n_val, m_val, V_val, c_val, p_val))
    end

    function _print_phase(label, phase_key, indices)
        isempty(indices) && return
        n_ph = _fmt(ustrip(us"mol",  state.n_phases[][phase_key]))
        m_ph = _fmt(ustrip(us"g",    state.m_phases[][phase_key]))
        V_ph = show_volume ?
               _fmt(ustrip(us"cm^3", state.V_phases[][phase_key])) : ""
        # Column sub-header — only show relevant columns per phase
        c_hdr = phase_key == :liquid ? "c [mol/L]" : ""
        p_hdr = phase_key == :gas    ? "p [bar]"   : ""
        println(io, _hhline())
        println(io, _row("# $label #", "n [mol]", "m [g]",
                         show_volume ? "V [cm³]" : "",
                         c_hdr, p_hdr))
        println(io, _htline())
        println(io, _row("tot. $label", n_ph, m_ph, V_ph, "", ""))
        println(io, _htline())
        for i in sort(indices; by=j->state.n[j], rev=true)
            _species_row(cs.species[i], state.n[i], phase_key)
        end
    end

    _print_phase("liquid", :liquid, cs.idx_aqueous)
    _print_phase("solid",  :solid,  cs.idx_crystal)
    _print_phase("gas",    :gas,    cs.idx_gas)

    # ── Total ─────────────────────────────────────────────────────────────────
    println(io, _hhline())
    println(io, _row("# TOTAL #", "n [mol]", "m [g]",
                     show_volume ? "V [cm³]" : "", "", ""))
    println(io, _htline())
    println(io, _row(
        "",
        _fmt(ustrip(us"mol",  state.n_phases[].total)),
        _fmt(ustrip(us"g",    state.m_phases[].total)),
        show_volume ? _fmt(ustrip(us"cm^3", state.V_phases[].total)) : "",
        "", "",
    ))
    println(io, _hhline())

    # ── Diagnostics ───────────────────────────────────────────────────────────
    any_diag = !isnothing(state.pH[])       ||
               !isnothing(state.pOH[])      ||
               !isnothing(state.porosity[]) ||
               !isnothing(state.saturation[])
    if any_diag
        isnothing(state.pH[])         || println(io, _full_row("pH",         _fmt4(state.pH[])))
        isnothing(state.pOH[])        || println(io, _full_row("pOH",        _fmt4(state.pOH[])))
        isnothing(state.porosity[])   || println(io, _full_row("porosity",   _fmt6(state.porosity[])))
        isnothing(state.saturation[]) || println(io, _full_row("saturation", _fmt6(state.saturation[])))
    end

    println(io, _botline())
end

"""
    Base.show(io::IO, state::ChemicalState)

Compact single-line representation of a `ChemicalState`.
"""
function Base.show(io::IO, state::ChemicalState)
    return print(io, "ChemicalState(T=", state.T[], ", P=", state.P[],
                 ", ", length(state.n), " species)")
end
