"""
    mutable struct ChemicalSystem{T<:AbstractSpecies, C, S} <: AbstractVector{T}

A mutable collection of chemical species with derived index structures
and stoichiometric matrices.

# Fields

  - `species`: ordered list of all species.
  - `dict_species`: fast O(1) lookup by species symbol.
  - `idx_aqueous`, `idx_crystal`, `idx_gas`: indices by aggregate state.
  - `idx_solutes`, `idx_solvent`, `idx_components`, `idx_gasfluid`: indices by class.
  - `CSM`: canonical stoichiometric matrix.
  - `SM`: stoichiometric matrix with respect to primaries.
"""
mutable struct ChemicalSystem{T<:AbstractSpecies, C, S} <: AbstractVector{T}
    species::Vector{T}
    dict_species::OrderedDict{String, T}

    # Indices by aggregate_state
    idx_aqueous::Vector{Int}
    idx_crystal::Vector{Int}
    idx_gas::Vector{Int}

    # Indices by class
    idx_solutes::Vector{Int}
    idx_solvent::Vector{Int}
    idx_components::Vector{Int}
    idx_gasfluid::Vector{Int}

    CSM::C
    SM::S
end

# ── Internal helpers ──────────────────────────────────────────────────────────

"""
    _reindex!(cs::ChemicalSystem) -> ChemicalSystem

Recompute all index vectors from the current species list.
Called internally after any mutation of `cs.species`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O"; aggregate_state=AS_AQUEOUS),
           Species("NaCl"; aggregate_state=AS_CRYSTAL),
       ]);

julia> cs.idx_aqueous
1-element Vector{Int64}:
 1

julia> cs.idx_crystal
1-element Vector{Int64}:
 2
```
"""
function _reindex!(cs::ChemicalSystem)
    idx(f) = findall(f, cs.species)
    cs.idx_aqueous    = idx(s -> aggregate_state(s) == AS_AQUEOUS)
    cs.idx_crystal    = idx(s -> aggregate_state(s) == AS_CRYSTAL)
    cs.idx_gas        = idx(s -> aggregate_state(s) == AS_GAS)
    cs.idx_solutes    = idx(s -> class(s) == SC_AQSOLUTE)
    cs.idx_solvent    = idx(s -> class(s) == SC_AQSOLVENT)
    cs.idx_components = idx(s -> class(s) == SC_COMPONENT)
    cs.idx_gasfluid   = idx(s -> class(s) == SC_GASFLUID)
    return cs
end

"""
    _rebuild_dict!(cs::ChemicalSystem) -> ChemicalSystem

Rebuild the symbol → species dictionary from scratch.
Called internally after any mutation of `cs.species`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS)]);

julia> haskey(cs.dict_species, "H2O")
true
```
"""
function _rebuild_dict!(cs::ChemicalSystem)
    empty!(cs.dict_species)
    for s in cs.species
        cs.dict_species[symbol(s)] = s
    end
    return cs
end

"""
    _rebuild_matrices!(cs::ChemicalSystem, primaries=cs.species) -> ChemicalSystem

Rebuild the canonical stoichiometric matrix `CSM` and the stoichiometric
matrix `SM` from the current species list.

Primaries default to all species, matching the behaviour of the main constructor.
Pass a custom primaries vector to preserve a specific choice across mutations.

Called internally after any mutation of `cs.species`.
"""
function _rebuild_matrices!(cs::ChemicalSystem, primaries=cs.species)
    cs.CSM = CanonicalStoichMatrix(cs.species)
    cs.SM  = StoichMatrix(cs.species, primaries)
    return cs
end

"""
    _update!(cs::ChemicalSystem, primaries=cs.species) -> ChemicalSystem

Full update of all derived fields after a mutation of `cs.species`:
  1. Rebuild the symbol → species dictionary.
  2. Recompute all aggregate-state and class index vectors.
  3. Rebuild `CSM` and `SM`.

All mutation methods (`setindex!`, `push!`, `append!`, `deleteat!`, `pop!`)
delegate to this function to keep the struct consistent.

!!! note
    Primaries default to all species. If custom primaries must survive
    mutations, store them as an extra field and pass them explicitly here.
"""
function _update!(cs::ChemicalSystem, primaries=cs.species)
    _rebuild_dict!(cs)       # keep dict_species consistent with species
    _reindex!(cs)            # keep index vectors consistent with species
    _rebuild_matrices!(cs, primaries)  # keep CSM and SM consistent with species
    return cs
end

# ── Constructors ──────────────────────────────────────────────────────────────

"""
    ChemicalSystem(species, primaries=species) -> ChemicalSystem

Construct a `ChemicalSystem` from a vector of species and an optional vector
of primary species used to build the stoichiometric matrix `SM`.

All derived fields (`dict_species`, index vectors, `CSM`, `SM`) are computed
at construction time.

# Arguments

  - `species`: vector of `AbstractSpecies`.
  - `primaries`: subset used as independent components (default: all species).

# Examples
```jldoctest
julia> sp = [
           Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("Na+"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
       ];

julia> cs = ChemicalSystem(sp);

julia> length(cs)
2

julia> cs["H2O"] == sp[1]
true
```
"""
function ChemicalSystem(
    species::AbstractVector{T},
    primaries::AbstractVector{<:AbstractSpecies}=species,
) where {T<:AbstractSpecies}
    idx(f) = findall(f, species)
    CSM = CanonicalStoichMatrix(species)
    SM  = StoichMatrix(species, primaries)
    return ChemicalSystem{T, typeof(CSM), typeof(SM)}(
        collect(species),                               # ensure owned Vector{T}
        OrderedDict(symbol(s) => s for s in species),  # symbol → species map
        idx(s -> aggregate_state(s) == AS_AQUEOUS),
        idx(s -> aggregate_state(s) == AS_CRYSTAL),
        idx(s -> aggregate_state(s) == AS_GAS),
        idx(s -> class(s) == SC_AQSOLUTE),
        idx(s -> class(s) == SC_AQSOLVENT),
        idx(s -> class(s) == SC_COMPONENT),
        idx(s -> class(s) == SC_GASFLUID),
        CSM,
        SM,
    )
end

"""
    ChemicalSystem(species, primaries::AbstractVector{<:AbstractString}) -> ChemicalSystem

Convenience constructor that resolves primary species from their symbol strings.

# Examples
```jldoctest
julia> sp = [
           Species("H2O";  aggregate_state=AS_AQUEOUS),
           Species("NaCl"; aggregate_state=AS_CRYSTAL),
       ];

julia> cs = ChemicalSystem(sp, ["H2O"]);

julia> length(cs)
2

julia> symbol.(cs.SM.primaries)
1-element Vector{String}:
 "H2O"
```
"""
function ChemicalSystem(
    species::AbstractVector{T},
    primaries::AbstractVector{<:AbstractString},
) where {T<:AbstractSpecies}
    # Resolve string symbols to species, preserving order
    primaries_species = @view species[symbol.(species) .∈ Ref(primaries)]
    return ChemicalSystem(species, primaries_species)
end

# ── AbstractVector interface ──────────────────────────────────────────────────

"""
    Base.size(cs::ChemicalSystem) -> Tuple

Return the size of the underlying species vector.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS)]);

julia> size(cs)
(1,)
```
"""
Base.size(cs::ChemicalSystem) = size(cs.species)

"""
    Base.getindex(cs::ChemicalSystem, i::Int) -> AbstractSpecies

Return the species at position `i`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS)]);

julia> cs[1] == Species("H2O"; aggregate_state=AS_AQUEOUS)
true
```
"""
Base.getindex(cs::ChemicalSystem, i::Int) = cs.species[i]

"""
    Base.getindex(cs::ChemicalSystem, i::AbstractString) -> AbstractSpecies

Return the species whose symbol matches `i`. Runs in O(1) via `dict_species`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS)]);

julia> cs["H2O"] == Species("H2O"; aggregate_state=AS_AQUEOUS)
true
```
"""
Base.getindex(cs::ChemicalSystem, i::AbstractString) = cs.dict_species[i]

"""
    Base.setindex!(cs::ChemicalSystem, x::AbstractSpecies, i::Int) -> ChemicalSystem

Replace the species at position `i` and update all derived fields,
including `CSM` and `SM` (rebuilt with `primaries = species`).

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS)]);

julia> cs[1] = Species("NaCl"; aggregate_state=AS_CRYSTAL);

julia> aggregate_state(cs[1]) == AS_CRYSTAL
true

julia> cs.idx_crystal
1-element Vector{Int64}:
 1

julia> isempty(cs.idx_aqueous)
true
```
"""
function Base.setindex!(cs::ChemicalSystem, x::AbstractSpecies, i::Int)
    cs.species[i] = x
    _update!(cs)  # rebuild dict, indices, CSM and SM
    return cs
end

# ── Mutation methods ──────────────────────────────────────────────────────────

"""
    Base.push!(cs::ChemicalSystem{T}, s::T) -> ChemicalSystem

Append a single species and update all derived fields,
including `CSM` and `SM`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS)]);

julia> push!(cs, Species("NaCl"; aggregate_state=AS_CRYSTAL));

julia> length(cs)
2

julia> cs.idx_crystal
1-element Vector{Int64}:
 2
```
"""
function Base.push!(cs::ChemicalSystem{T}, s::T) where {T<:AbstractSpecies}
    push!(cs.species, s)
    _update!(cs)  # rebuild dict, indices, CSM and SM
    return cs
end

"""
    Base.append!(cs::ChemicalSystem{T}, v::AbstractVector{<:AbstractSpecies}) -> ChemicalSystem

Append multiple species at once and update all derived fields.
More efficient than successive `push!` calls since `_update!` is called only once,
rebuilding `CSM` and `SM` a single time.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS)]);

julia> append!(cs, [
           Species("NaCl"; aggregate_state=AS_CRYSTAL),
           Species("CO2";  aggregate_state=AS_GAS),
       ]);

julia> length(cs)
3

julia> length(cs.idx_aqueous), length(cs.idx_crystal), length(cs.idx_gas)
(1, 1, 1)
```
"""
function Base.append!(
    cs::ChemicalSystem{T},
    v::AbstractVector{<:AbstractSpecies},
) where {T<:AbstractSpecies}
    append!(cs.species, v)
    _update!(cs)  # single rebuild for all appended species
    return cs
end

"""
    Base.deleteat!(cs::ChemicalSystem, i) -> ChemicalSystem

Remove the species at index `i` and update all derived fields,
including `CSM` and `SM`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O";  aggregate_state=AS_AQUEOUS),
           Species("NaCl"; aggregate_state=AS_CRYSTAL),
       ]);

julia> deleteat!(cs, 2);

julia> length(cs)
1

julia> isempty(cs.idx_crystal)
true
```
"""
function Base.deleteat!(cs::ChemicalSystem, i)
    deleteat!(cs.species, i)
    _update!(cs)  # rebuild dict, indices, CSM and SM after removal
    return cs
end

"""
    Base.pop!(cs::ChemicalSystem) -> AbstractSpecies

Remove and return the last species, updating all derived fields,
including `CSM` and `SM`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O";  aggregate_state=AS_AQUEOUS),
           Species("NaCl"; aggregate_state=AS_CRYSTAL),
       ]);

julia> s = pop!(cs);

julia> symbol(s)
"NaCl"

julia> length(cs)
1

julia> isempty(cs.idx_crystal)
true
```
"""
function Base.pop!(cs::ChemicalSystem)
    s = pop!(cs.species)     # remove last species and capture it
    _update!(cs)             # rebuild dict, indices, CSM and SM
    return s                 # return the removed species to the caller
end

# ── Merge ─────────────────────────────────────────────────────────────────────

"""
    Base.merge(cs1::ChemicalSystem, cs2::ChemicalSystem) -> ChemicalSystem

Merge two `ChemicalSystem` instances into a single one.

Species are unioned by symbol — duplicates from `cs2` are discarded.
`CSM` and `SM` are rebuilt from scratch from the full species list.
Primaries are taken as the union of both systems' primaries, filtered
to those actually present in the merged species list.

In case of symbol conflict, `cs1` takes priority over `cs2`.

# Examples
```jldoctest
julia> cs1 = ChemicalSystem([
           Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("H+";  aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
       ]);

julia> cs2 = ChemicalSystem([
           Species("OH-"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
           Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
       ]);

julia> cs = merge(cs1, cs2);

julia> length(cs)
3

julia> symbol.(cs.SM.primaries)
2-element Vector{String}:
 "H2O"
 "H+"
```
"""
function Base.merge(cs1::ChemicalSystem, cs2::ChemicalSystem)
    # Build a set of existing symbols from cs1 for fast duplicate detection
    existing_symbols = Set(symbol.(cs1.species))

    # Append only species from cs2 whose symbol is not already in cs1
    extra_species = filter(s -> symbol(s) ∉ existing_symbols, cs2.species)
    all_species   = vcat(cs1.species, extra_species)

    # Union of primaries: cs1 first, then new ones from cs2 not already present
    existing_primary_symbols = Set(symbol.(cs1.SM.primaries))
    extra_primaries = filter(
        p -> symbol(p) ∉ existing_primary_symbols,
        cs2.SM.primaries,
    )
    all_primaries = vcat(cs1.SM.primaries, extra_primaries)

    # Drop primaries that are not present in the merged species list
    # (can happen if a cs2 primary was shadowed or not included)
    all_symbols   = Set(symbol.(all_species))
    all_primaries = filter(p -> symbol(p) ∈ all_symbols, all_primaries)

    # Delegate to main constructor — rebuilds CSM and SM from scratch
    return ChemicalSystem(all_species, all_primaries)
end

"""
    Base.merge(css::ChemicalSystem...) -> ChemicalSystem

Merge an arbitrary number of `ChemicalSystem` instances left-to-right.
Earlier systems take priority over later ones in case of symbol conflicts.

# Examples
```jldoctest
julia> cs1 = ChemicalSystem([Species("H2O";  aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> cs2 = ChemicalSystem([Species("H+";   aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE)]);

julia> cs3 = ChemicalSystem([Species("OH-";  aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE)]);

julia> cs = merge(cs1, cs2, cs3);

julia> length(cs)
3
```
"""
Base.merge(css::ChemicalSystem...) = reduce(merge, css)

# ── Views by aggregate state ──────────────────────────────────────────────────

"""
    aqueous(cs::ChemicalSystem) -> SubArray

Return a view of all aqueous species.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O";  aggregate_state=AS_AQUEOUS),
           Species("NaCl"; aggregate_state=AS_CRYSTAL),
       ]);

julia> length(aqueous(cs))
1

julia> aggregate_state(aqueous(cs)[1]) == AS_AQUEOUS
true
```
"""
aqueous(cs::ChemicalSystem) = @view cs.species[cs.idx_aqueous]

"""
    crystal(cs::ChemicalSystem) -> SubArray

Return a view of all crystalline species.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O";  aggregate_state=AS_AQUEOUS),
           Species("NaCl"; aggregate_state=AS_CRYSTAL),
       ]);

julia> aggregate_state(crystal(cs)[1]) == AS_CRYSTAL
true
```
"""
crystal(cs::ChemicalSystem) = @view cs.species[cs.idx_crystal]

"""
    gas(cs::ChemicalSystem) -> SubArray

Return a view of all gas-phase species.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("CO2"; aggregate_state=AS_GAS)]);

julia> aggregate_state(gas(cs)[1]) == AS_GAS
true
```
"""
gas(cs::ChemicalSystem) = @view cs.species[cs.idx_gas]

# ── Views by class ────────────────────────────────────────────────────────────

"""
    solutes(cs::ChemicalSystem) -> SubArray

Return a view of all aqueous solute species.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("Na+"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE)]);

julia> class(solutes(cs)[1]) == SC_AQSOLUTE
true
```
"""
solutes(cs::ChemicalSystem) = @view cs.species[cs.idx_solutes]

"""
    solvent(cs::ChemicalSystem) -> AbstractSpecies

Return the unique solvent species directly (not a view),
since a chemical system contains at most one solvent.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> class(solvent(cs)) == SC_AQSOLVENT
true
```
"""
solvent(cs::ChemicalSystem) = cs.species[cs.idx_solvent][1]  # unique element, return directly

"""
    components(cs::ChemicalSystem) -> SubArray

Return a view of all component species.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("SiO2"; aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)]);

julia> class(components(cs)[1]) == SC_COMPONENT
true
```
"""
components(cs::ChemicalSystem) = @view cs.species[cs.idx_components]

"""
    gasfluid(cs::ChemicalSystem) -> SubArray

Return a view of all gas/fluid species.

# Examples
```jldoctest
julia> cs = ChemicalSystem([Species("CO2"; aggregate_state=AS_GAS, class=SC_GASFLUID)]);

julia> class(gasfluid(cs)[1]) == SC_GASFLUID
true
```
"""
gasfluid(cs::ChemicalSystem) = @view cs.species[cs.idx_gasfluid]
