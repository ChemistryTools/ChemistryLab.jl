using DynamicQuantities

# ── Abstract surface area model ───────────────────────────────────────────────

"""
    abstract type AbstractSurfaceModel end

Base type for models that compute the reactive surface area of a mineral phase.

Concrete subtypes must implement:
```julia
surface_area(model, n::Real, molar_mass::Real) -> Real
```
returning the reactive surface area in m².

All methods must be AD-compatible (no `Float64` casts).
"""
abstract type AbstractSurfaceModel end

# ── FixedSurfaceArea ──────────────────────────────────────────────────────────

"""
    struct FixedSurfaceArea{T<:Real} <: AbstractSurfaceModel

Constant reactive surface area, independent of mineral abundance.

Suitable for short simulations or when the surface area is externally controlled
(e.g. from BET measurements on a fixed mass of powder).

# Fields

  - `A`: total reactive surface area [m²].

# Examples

```julia
FixedSurfaceArea(0.5)    # 0.5 m²
```
"""
struct FixedSurfaceArea{T <: Real} <: AbstractSurfaceModel
    A::T
end

"""
    surface_area(model::FixedSurfaceArea, n::Real, molar_mass::Real) -> Real

Return the fixed surface area `model.A` [m²], independent of moles `n`.
AD-compatible.
"""
surface_area(model::FixedSurfaceArea, ::Real, ::Real) = model.A

# ── BETSurfaceArea ────────────────────────────────────────────────────────────

"""
    struct BETSurfaceArea{T<:Real} <: AbstractSurfaceModel

Reactive surface area that scales with the mineral mass, following a BET
(Brunauer-Emmett-Teller) specific-surface-area measurement.

```
A = A_spec × n × M_mineral       [m²]
```

where `A_spec` [m²/kg] is the specific BET surface area, `n` is the current
molar amount [mol], and `M_mineral` is the molar mass [kg/mol].

This is the standard approach in reactive-transport models (Palandri & Kharaka 2004).

# Fields

  - `A_specific`: specific BET surface area [m²/kg].

# Examples

```julia
# Calcite: 0.09 m²/g = 90 m²/kg
BETSurfaceArea(90.0)
```
"""
struct BETSurfaceArea{T <: Real} <: AbstractSurfaceModel
    A_specific::T   # m²/kg
end

"""
    surface_area(model::BETSurfaceArea, n::Real, molar_mass::Real) -> Real

Return `A_specific × n × molar_mass` [m²].  Clamps at zero to avoid negative
surface areas when `n → 0`. AD-compatible.
"""
function surface_area(model::BETSurfaceArea, n::Real, molar_mass::Real)
    return model.A_specific * max(n, zero(n)) * molar_mass
end

# ── KineticReaction ───────────────────────────────────────────────────────────

"""
    struct KineticReaction{R, M<:AbstractRateModel, S<:AbstractSurfaceModel, H}

Associates a mineral species (or thermodynamic [`Reaction`](@ref)) with a
kinetic [`AbstractRateModel`](@ref) and an [`AbstractSurfaceModel`](@ref).

Kinetic reactions are "slow" reactions whose extent is governed by a rate law
rather than by instantaneous equilibrium. Typically these are mineral
dissolution/precipitation reactions.

# Fields

  - `reaction`: the underlying [`Reaction`](@ref) / [`CemReaction`](@ref) *or* a
    bare [`Species`](@ref) (when constructed via the convenience constructor).
    Used only to retrieve the molar mass for [`BETSurfaceArea`](@ref) calculations.
  - `rate_model`: an [`AbstractRateModel`](@ref) computing r [mol/s].
  - `surface_model`: an [`AbstractSurfaceModel`](@ref) computing A [m²].
  - `idx_mineral`: index of the mineral species in the parent `ChemicalSystem`.
  - `stoich`: stoichiometric coefficient vector for all species in the system.
    Sign convention: positive for products, negative for reactants.
    Used to evaluate the saturation ratio Ω (ignored by
    [`ParrotKillohRateModel`](@ref)).
  - `heat_per_mol`: enthalpy of reaction [J/mol], positive = exothermic (heat released).
    Used by calorimeter models to compute q̇(t). When `nothing` (default), the
    enthalpy is derived from the stoichiometric sum of species `:ΔₐH⁰` values —
    only meaningful when `reaction isa AbstractReaction` with full stoichiometry.
    For convenience-constructed reactions (bare [`Species`](@ref)) without explicit
    reaction stoichiometry, set this to a literature value (e.g. 114_634.0 J/mol
    for C₃S hydration).

# Constructors

The **recommended** user-facing constructor takes a [`ChemicalSystem`](@ref)
and a species name — no manual index lookup or stoichiometry needed:

```julia
# Convenience: system + name → index and stoichiometry built automatically
kr = KineticReaction(cs, "C3S", PK_PARAMS_C3S, FixedSurfaceArea(1.0))

# With explicit heat of hydration for calorimetry [J/mol]
kr = KineticReaction(cs, "C3S", PK_PARAMS_C3S, FixedSurfaceArea(1.0);
                     heat_per_mol = 114_634.0)
```

The **low-level** constructor accepts an explicit `Reaction` (or `Species`),
index, and stoichiometry vector for cases where full dissolution-reaction
stoichiometry must be supplied (e.g. [`TransitionStateRateModel`](@ref)):

```julia
k_acid    = arrhenius_rate_constant(5.012e-1, 14400.0)
k_neutral = arrhenius_rate_constant(1.549e-6, 23500.0)

model = TransitionStateRateModel([
    RateMechanism(k_acid,    1.0, 1.0, [RateModelCatalyst("H+", 1.0)]),
    RateMechanism(k_neutral, 1.0, 1.0),
])

kr = KineticReaction(
    reaction_calcite,           # Reaction object from ChemicalSystem
    model,
    BETSurfaceArea(90.0),       # 90 m²/kg (BET)
    idx_calcite,                # integer index in ChemicalSystem
    stoich_calcite,             # stoichiometric row vector
)
```
"""
struct KineticReaction{R, M <: AbstractRateModel, S <: AbstractSurfaceModel, H}
    reaction::R        # AbstractReaction or AbstractSpecies (for convenience ctor)
    rate_model::M
    surface_model::S
    idx_mineral::Int
    stoich::Vector{Float64}    # stoich coefficients for all species in system
    heat_per_mol::H            # Nothing or Float64: enthalpy of reaction [J/mol], positive = exothermic

    # Inner constructor: always validates, regardless of which outer constructor is used.
    function KineticReaction{R, M, S, H}(
            reaction::R,
            rate_model::M,
            surface_model::S,
            idx_mineral::Int,
            stoich::Vector{Float64},
            heat_per_mol::H,
        ) where {R, M <: AbstractRateModel, S <: AbstractSurfaceModel, H}
        idx_mineral > 0 || throw(ArgumentError("idx_mineral must be a positive integer"))
        isempty(stoich) && throw(ArgumentError("stoich cannot be empty"))
        return new{R, M, S, H}(reaction, rate_model, surface_model, idx_mineral, stoich, heat_per_mol)
    end
end

"""
    KineticReaction(reaction, rate_model, surface_model, idx_mineral, stoich;
                    heat_per_mol=nothing)

Low-level constructor: explicit `Reaction` (or `Species`), index, and stoichiometry.
Validates that `idx_mineral > 0` and `stoich` is non-empty.

The first argument accepts any type (typically `AbstractReaction` or `AbstractSpecies`);
the correct [`molar_mass`](@ref) dispatch is resolved at call time based on the stored type.

`heat_per_mol` [J/mol] overrides stoichiometric enthalpy calculation for calorimetry.
"""
function KineticReaction(
        reaction::R,
        rate_model::M,
        surface_model::S,
        idx_mineral::Integer,
        stoich::AbstractVector{<:Real};
        heat_per_mol::Union{Nothing, Float64} = nothing,
    ) where {R, M <: AbstractRateModel, S <: AbstractSurfaceModel}
    return KineticReaction{R, M, S, typeof(heat_per_mol)}(
        reaction, rate_model, surface_model, Int(idx_mineral), Float64.(stoich), heat_per_mol
    )
end

"""
    KineticReaction(cs::ChemicalSystem, species_name, rate_model, surface_model;
                    stoich=nothing) -> KineticReaction

Convenience constructor: look up `species_name` in `cs`, determine the mineral
index and default stoichiometry automatically.

The default stoichiometry is a vector of zeros with `-1.0` at the mineral index,
meaning "pure mineral dissolution" (no products tracked). This is correct for
[`ParrotKillohRateModel`](@ref) (which ignores Ω) and for any model where you
only care about the mineral mole balance, not the full reaction stoichiometry.

For [`TransitionStateRateModel`](@ref) with explicit catalyst/product effects,
supply the full stoichiometry via the `stoich` keyword.

# Arguments

  - `cs`: the [`ChemicalSystem`](@ref) containing the mineral.
  - `species_name`: PHREEQC formula or symbol string (e.g. `"C3S"`, `"Calcite"`).
  - `rate_model`: an [`AbstractRateModel`](@ref).
  - `surface_model`: an [`AbstractSurfaceModel`](@ref).
  - `stoich`: optional stoichiometric coefficient vector (length = `length(cs.species)`).

# Examples

```julia
using ChemistryLab, OrdinaryDiffEq

# Build kinetic reactions for all four clinker phases — no manual indexing
kr_C3S  = KineticReaction(cs, "C3S",  PK_PARAMS_C3S,  FixedSurfaceArea(1.0))
kr_C2S  = KineticReaction(cs, "C2S",  PK_PARAMS_C2S,  FixedSurfaceArea(1.0))
kr_C3A  = KineticReaction(cs, "C3A",  PK_PARAMS_C3A,  FixedSurfaceArea(1.0))
kr_C4AF = KineticReaction(cs, "C4AF", PK_PARAMS_C4AF, FixedSurfaceArea(1.0))

kp  = KineticsProblem(cs, [kr_C3S, kr_C2S, kr_C3A, kr_C4AF], state0, (0.0, 7*86400.0))
sol = integrate(kp)
```
"""
function KineticReaction(
        cs::ChemicalSystem,
        species_name::AbstractString,
        rate_model::AbstractRateModel,
        surface_model::AbstractSurfaceModel;
        stoich::Union{Nothing, AbstractVector{<:Real}} = nothing,
        heat_per_mol::Union{Nothing, Float64} = nothing,
    )
    # Find species index
    idx = findfirst(
        sp -> phreeqc(formula(sp)) == species_name || string(symbol(sp)) == species_name,
        cs.species,
    )
    isnothing(idx) && throw(
        ArgumentError(
            "Species \"$species_name\" not found in ChemicalSystem. " *
                "Use phreeqc(formula(sp)) or symbol(sp) to check species names."
        )
    )

    # Default stoichiometry: -1 for the mineral, 0 for all others
    stoich_vec = if isnothing(stoich)
        s = zeros(Float64, length(cs.species))
        s[idx] = -1.0
        s
    else
        length(stoich) == length(cs.species) || throw(
            DimensionMismatch(
                "stoich length $(length(stoich)) ≠ number of species $(length(cs.species))"
            )
        )
        Float64.(stoich)
    end

    sp = cs.species[idx]
    return KineticReaction{typeof(sp), typeof(rate_model), typeof(surface_model), typeof(heat_per_mol)}(
        sp, rate_model, surface_model, Int(idx), stoich_vec, heat_per_mol
    )
end

"""
    molar_mass(kr::KineticReaction) -> Float64

Return the molar mass of the mineral species [kg/mol], used internally for
[`BETSurfaceArea`](@ref) calculations.

Two dispatch paths:
- If `kr.reaction isa AbstractReaction`: searches `reaction.reactants` for a
  species with an `:M` property.
- If `kr.reaction isa AbstractSpecies`: reads `:M` directly from the species
  (set by the convenience constructor [`KineticReaction`](@ref)).

Falls back to `0.1` kg/mol when `:M` is unavailable.
"""
function molar_mass(kr::KineticReaction{<:AbstractReaction})
    # Search reactants for a species that carries a molar mass property
    for (sp_obj, _) in kr.reaction.reactants
        if haskey(properties(sp_obj), :M)
            return ustrip(us"kg/mol", sp_obj[:M])
        end
    end
    return 0.1   # fallback: ~100 g/mol
end

function molar_mass(kr::KineticReaction{<:AbstractSpecies})
    sp = kr.reaction   # the bare Species stored by the convenience constructor
    haskey(properties(sp), :M) && return ustrip(us"kg/mol", sp[:M])
    return 0.1
end
