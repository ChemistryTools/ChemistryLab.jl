using LinearAlgebra
using OrderedCollections

"""
    struct ChemicalSystem{T<:AbstractSpecies, R<:AbstractReaction, C, S, SS} <: AbstractVector{T}

An immutable, fully typed collection of chemical species and reactions
with derived index structures and stoichiometric matrices.

Immutability guarantees that all derived fields (`dict_species`, `dict_reactions`,
index vectors, `CSM`, `SM`) remain consistent with `species` and `reactions`
throughout the lifetime of the object. To modify the system, use `merge` to
construct a new `ChemicalSystem`.

# Fields

  - `species`: ordered list of all species.
  - `dict_species`: fast O(1) lookup by species symbol.
  - `idx_aqueous`, `idx_crystal`, `idx_gas`: indices by aggregate state.
  - `idx_solutes`, `idx_solvent`, `idx_components`, `idx_gasfluid`: indices by class.
  - `reactions`: ordered list of all reactions.
  - `dict_reactions`: fast O(1) lookup by reaction symbol.
  - `CSM`: canonical stoichiometric matrix.
  - `SM`: stoichiometric matrix with respect to primaries.
  - `solid_solutions`: `Nothing` when no solid solutions are present, or a concrete
    `Vector{<:AbstractSolidSolutionPhase}` describing each solid-solution phase and
    its end-members. Populated via the `solid_solutions` keyword constructor.
  - `ss_groups`: for each solid solution, the indices of its end-members in `species`.
  - `idx_ssendmembers`: union of all end-member indices (flattened `ss_groups`).
"""
struct ChemicalSystem{T <: AbstractSpecies, R <: AbstractReaction, C, S, SS} <:
    AbstractVector{T}
    species::Vector{T}
    dict_species::Dict{String, T}               # fast O(1) lookup by symbol

    # Indices by aggregate_state
    idx_aqueous::Vector{Int}
    idx_crystal::Vector{Int}
    idx_gas::Vector{Int}

    # Indices by class
    idx_solutes::Vector{Int}
    idx_solvent::Vector{Int}
    idx_components::Vector{Int}
    idx_gasfluid::Vector{Int}

    reactions::Vector{R}
    dict_reactions::Dict{String, R}             # fast O(1) lookup by reaction symbol

    CSM::C                                      # canonical stoichiometric matrix — typed for performance
    SM::S                                       # stoichiometric matrix w.r.t. primaries — typed for performance

    # Solid solutions — SS = Nothing (no SS) or Vector{<:AbstractSolidSolutionPhase}
    solid_solutions::SS
    ss_groups::Vector{Vector{Int}}              # per-SS end-member indices
    idx_ssendmembers::Vector{Int}               # all end-member indices (flattened)
end

# ── Constructors ──────────────────────────────────────────────────────────────

"""
    ChemicalSystem(species, primaries=species; reactions, solid_solutions=nothing) -> ChemicalSystem

Construct a fully typed `ChemicalSystem` from a vector of species,
an optional vector of primary species, optional reactions, and optional solid-solution phases.

All derived fields are computed once at construction time and remain
consistent for the lifetime of the object.

# Arguments

  - `species`: vector of `AbstractSpecies`.
  - `primaries`: subset used as independent components (default: all species).
  - `reactions`: vector of `AbstractReaction` (default: empty).
  - `solid_solutions`: vector of [`SolidSolutionPhase`](@ref) (default: `nothing`).
    When provided, end-members must already appear in `species` (matched by symbol) and
    must carry `aggregate_state = AS_CRYSTAL` and `class = SC_SSENDMEMBER`.

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

```jldoctest
julia> em1 = Species("AFm1"; aggregate_state=AS_CRYSTAL, class=SC_SSENDMEMBER);

julia> em2 = Species("AFm2"; aggregate_state=AS_CRYSTAL, class=SC_SSENDMEMBER);

julia> ss = SolidSolutionPhase("AFm", [em1, em2]);

julia> cs = ChemicalSystem([em1, em2]; solid_solutions=[ss]);

julia> cs.ss_groups
1-element Vector{Vector{Int64}}:
 [1, 2]

julia> cs.idx_ssendmembers
2-element Vector{Int64}:
 1
 2
```
"""
function ChemicalSystem(
        species::AbstractVector{T},
        primaries::AbstractVector{<:AbstractSpecies} = species;
        reactions::AbstractVector{R} = AbstractReaction[],
        solid_solutions::Union{Nothing, AbstractVector{<:AbstractSolidSolutionPhase}} = nothing,
    ) where {T <: AbstractSpecies, R <: AbstractReaction}
    idx(f) = findall(f, species)
    CSM = CanonicalStoichMatrix(species)
    SM = StoichMatrix(species, primaries)

    if isnothing(solid_solutions)
        return ChemicalSystem{T, R, typeof(CSM), typeof(SM), Nothing}(
            collect(T, species),                                    # owned Vector{T}
            Dict{String, T}(symbol(s) => s for s in species),      # symbol → species map
            idx(s -> aggregate_state(s) == AS_AQUEOUS),
            idx(s -> aggregate_state(s) == AS_CRYSTAL),
            idx(s -> aggregate_state(s) == AS_GAS),
            idx(s -> class(s) == SC_AQSOLUTE),
            idx(s -> class(s) == SC_AQSOLVENT),
            idx(s -> class(s) == SC_COMPONENT),
            idx(s -> class(s) == SC_GASFLUID),
            collect(R, reactions),                                  # owned Vector{R}
            Dict{String, R}(symbol(r) => r for r in reactions),    # symbol → reaction map
            CSM,
            SM,
            nothing,        # solid_solutions — SS = Nothing
            Vector{Int}[],  # ss_groups
            Int[],          # idx_ssendmembers
        )
    else
        # Build ss_groups: for each solid solution, find the indices of its end-members
        # in the species vector (by symbol matching).
        ss_groups = map(solid_solutions) do ss
            map(end_members(ss)) do em
                idx_em = findfirst(s -> symbol(s) == symbol(em), species)
                idx_em === nothing &&
                    error(
                    "SolidSolutionPhase \"$(name(ss))\": end-member \"$(symbol(em))\" " *
                        "not found in the species list. Add it to the species vector first.",
                )
                idx_em
            end
        end
        idx_ssendmembers = isempty(ss_groups) ? Int[] : vcat(ss_groups...)
        ss = collect(solid_solutions)

        return ChemicalSystem{T, R, typeof(CSM), typeof(SM), typeof(ss)}(
            collect(T, species),                                    # owned Vector{T}
            Dict{String, T}(symbol(s) => s for s in species),      # symbol → species map
            idx(s -> aggregate_state(s) == AS_AQUEOUS),
            idx(s -> aggregate_state(s) == AS_CRYSTAL),
            idx(s -> aggregate_state(s) == AS_GAS),
            idx(s -> class(s) == SC_AQSOLUTE),
            idx(s -> class(s) == SC_AQSOLVENT),
            idx(s -> class(s) == SC_COMPONENT),
            idx(s -> class(s) == SC_GASFLUID),
            collect(R, reactions),                                  # owned Vector{R}
            Dict{String, R}(symbol(r) => r for r in reactions),    # symbol → reaction map
            CSM,
            SM,
            ss,
            ss_groups,
            idx_ssendmembers,
        )
    end
end

"""
    ChemicalSystem(species, primaries::AbstractVector{<:AbstractString}; reactions) -> ChemicalSystem

Convenience constructor that resolves primary species from their symbol strings.

# Examples
```jldoctest
julia> sp = [
           Species("H2O";  aggregate_state=AS_AQUEOUS),
           Species("NaCl"; aggregate_state=AS_CRYSTAL),
       ];

julia> cs = ChemicalSystem(sp, ["H2O"]);

julia> symbol.(cs.SM.primaries)
1-element Vector{String}:
 "H2O"
```
"""
function ChemicalSystem(
        species::AbstractVector{T},
        primaries::AbstractVector{<:AbstractString};
        reactions::AbstractVector{R} = AbstractReaction[],
        solid_solutions::Union{Nothing, AbstractVector{<:AbstractSolidSolutionPhase}} = nothing,
    ) where {T <: AbstractSpecies, R <: AbstractReaction}
    # Resolve string symbols to species objects, preserving order
    primaries_species = species[symbol.(species) .∈ Ref(primaries)]
    return ChemicalSystem(
        species, primaries_species; reactions = reactions, solid_solutions = solid_solutions
    )
end

# ── Solid solution accessor ───────────────────────────────────────────────────

"""
    solid_solutions(cs::ChemicalSystem) -> Nothing | Vector{<:AbstractSolidSolutionPhase}

Return the registered solid-solution phases, or `nothing` if none were declared.

# Examples
```jldoctest
julia> em1 = Species("Em1"; aggregate_state=AS_CRYSTAL, class=SC_SSENDMEMBER);

julia> em2 = Species("Em2"; aggregate_state=AS_CRYSTAL, class=SC_SSENDMEMBER);

julia> cs = ChemicalSystem(
           [em1, em2];
           solid_solutions=[SolidSolutionPhase("SS", [em1, em2])],
       );

julia> solid_solutions(cs) isa Vector
true

julia> length(solid_solutions(cs))
1
```
"""
solid_solutions(cs::ChemicalSystem) = cs.solid_solutions

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

# ── Reaction accessor ─────────────────────────────────────────────────────────

"""
    get_reaction(cs::ChemicalSystem, sym::AbstractString) -> AbstractReaction

Return the reaction identified by symbol `sym`. Runs in O(1) via `dict_reactions`.

# Examples
```jldoctest
julia> cs = ChemicalSystem(
           [Species("H2O"; aggregate_state=AS_AQUEOUS)];
           reactions=[Reaction("H2O = H+ + OH-"; symbol="water_diss")],
       );

julia> symbol(get_reaction(cs, "water_diss"))
"water_diss"
```
"""
get_reaction(cs::ChemicalSystem, sym::AbstractString) = cs.dict_reactions[sym]

# ── Merge ─────────────────────────────────────────────────────────────────────

"""
    Base.merge(cs1::ChemicalSystem, cs2::ChemicalSystem) -> ChemicalSystem

Construct a new `ChemicalSystem` from the union of two systems.

Species and reactions are unioned by symbol — duplicates from `cs2` are discarded.
`CSM` and `SM` are built from scratch from the full species list.
Primaries are taken as the union of both systems' primaries, filtered
to those actually present in the merged species list.

In case of symbol conflict (species or reactions), `cs1` takes priority over `cs2`.
The return type is inferred from the merged collections and may differ from
`typeof(cs1)` or `typeof(cs2)` if they contain different concrete types.

# Examples
```jldoctest
julia> cs1 = ChemicalSystem(
           [Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
            Species("H+";  aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE)];
           reactions=[Reaction("H2O = H+ + OH-"; symbol="water_diss")],
       );

julia> cs2 = ChemicalSystem(
           [Species("OH-"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
            Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)];
           reactions=[Reaction("CO2 + H2O = H2CO3"; symbol="co2_hyd")],
       );

julia> cs = merge(cs1, cs2);

julia> length(cs)
3

julia> length(cs.reactions)
2
```
"""
function Base.merge(cs1::ChemicalSystem, cs2::ChemicalSystem)
    # Build lookup sets for fast duplicate detection
    existing_symbols = Set(symbol.(cs1.species))

    # Append only species from cs2 not already present in cs1 (cs1 wins on conflict)
    extra_species = filter(s -> symbol(s) ∉ existing_symbols, cs2.species)
    all_species = vcat(cs1.species, extra_species)

    # Union of primaries: cs1 first, then new ones from cs2 not already present
    existing_primary_symbols = Set(symbol.(cs1.SM.primaries))
    extra_primaries = filter(
        p -> symbol(p) ∉ existing_primary_symbols,
        cs2.SM.primaries,
    )
    all_primaries = vcat(cs1.SM.primaries, extra_primaries)

    # Drop primaries absent from the merged species list
    all_species_symbols = Set(symbol.(all_species))
    all_primaries = filter(p -> symbol(p) ∈ all_species_symbols, all_primaries)

    # Union of reactions by symbol — cs1 wins on conflict
    existing_reaction_symbols = Set(symbol.(cs1.reactions))
    extra_reactions = filter(r -> symbol(r) ∉ existing_reaction_symbols, cs2.reactions)
    all_reactions = vcat(cs1.reactions, extra_reactions)

    # Construct a new ChemicalSystem — all derived fields rebuilt from scratch
    return ChemicalSystem(all_species, all_primaries; reactions = all_reactions)
end

"""
    Base.merge(css::ChemicalSystem...) -> ChemicalSystem

Construct a new `ChemicalSystem` from the union of an arbitrary number of systems,
processed left-to-right. Earlier systems take priority over later ones
in case of symbol conflicts.

# Examples
```jldoctest
julia> cs1 = ChemicalSystem([Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)]);

julia> cs2 = ChemicalSystem([Species("H+"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE)]);

julia> cs3 = ChemicalSystem([Species("OH-"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE)]);

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
