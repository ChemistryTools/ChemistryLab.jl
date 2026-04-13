using DynamicQuantities
using LinearAlgebra
using SciMLBase

# ── KineticsProblem ───────────────────────────────────────────────────────────

"""
    struct KineticsProblem{CS, KR, ES, AM<:AbstractActivityModel}

Encapsulates the full specification of a kinetics simulation:
chemical system, kinetic reactions, initial state, time span, and solver ingredients.

Construct via [`KineticsProblem`](@ref) (outer constructor).

# Fields

  - `system`: the [`ChemicalSystem`](@ref) holding all species and reactions.
  - `kinetic_reactions`: vector of [`KineticReaction`](@ref) objects.
  - `initial_state`: [`ChemicalState`](@ref) providing initial mole amounts, T, P.
  - `tspan`: `(t_start, t_end)` time interval [s].
  - `activity_model`: [`AbstractActivityModel`](@ref) used at every ODE step to
    compute log-activities.
  - `equilibrium_solver`: an [`EquilibriumSolver`](@ref) re-equilibrating the fast
    (aqueous) species at every ODE evaluation. `nothing` disables re-speciation.
  - `S_kin`: stoichiometric submatrix `[n_reactions × n_species]` — row `i` is the
    stoichiometric coefficient vector for kinetic reaction `i`.
  - `idx_kin`: integer indices of each kinetic reaction's controlling mineral in `system`
    (length = n_reactions, may contain duplicates when several reactions share a mineral).
  - `idx_kin_unique`: sorted unique kinetic species indices — defines the ODE state vector
    `u`. Includes the rate-controlling minerals (`idx_kin`) PLUS all AS_CRYSTAL species
    with non-zero stoich in any kinetic reaction (Leal et al. 2015 kinetic classification).

# Usage

```julia
kp = KineticsProblem(cs, [kr_calcite], state0, (0.0, 3600.0))
```

See also: [`integrate`](@ref), [`KineticsSolver`](@ref).
"""
struct KineticsProblem{
        CS <: ChemicalSystem,
        KR <: AbstractVector,
        ES,
        AM <: AbstractActivityModel,
    }
    system::CS
    kinetic_reactions::KR
    initial_state::ChemicalState
    tspan::Tuple{Float64, Float64}
    activity_model::AM
    equilibrium_solver::ES
    S_kin::Matrix{Float64}
    idx_kin::Vector{Int}
    idx_kin_unique::Vector{Int}
end

"""
    KineticsProblem(system, kinetic_reactions, initial_state, tspan;
                    activity_model = DiluteSolutionModel(),
                    equilibrium_solver = nothing) -> KineticsProblem

Construct a [`KineticsProblem`](@ref).

Multiple [`KineticReaction`](@ref) objects may share the same mineral (`idx_mineral`),
enabling multi-pathway kinetics. The ODE state has **one entry per unique mineral**;
contributions from all reactions are accumulated.

# Arguments

  - `tspan`: `(t0, tf)` time interval. Plain `Real` → SI [s]; `Quantity` → converted
    (e.g. `(0.0u"s", 7.0u"d")` or `(0.0, 168.0u"hr")`). Mixed tuples supported.

# Examples

```julia
pk = parrot_killoh(PK_PARAMS_C3S, "C3S")
kr = KineticReaction(cs, "C3S", pk; heat_per_mol = 114_634.0)
kp = KineticsProblem(cs, [kr], state0, (0.0, 7 * 86400.0))
```
"""
function KineticsProblem(
        system::ChemicalSystem,
        kinetic_reactions::AbstractVector,
        initial_state::ChemicalState,
        tspan::Tuple;
        activity_model::AbstractActivityModel = DiluteSolutionModel(),
        equilibrium_solver = nothing,
    )
    n_sp = length(system.species)
    n_rxn = length(kinetic_reactions)

    idx_kin = Int[kr.idx_mineral for kr in kinetic_reactions]
    S_kin = zeros(Float64, n_rxn, n_sp)
    for (i, kr) in enumerate(kinetic_reactions)
        length(kr.stoich) == n_sp || throw(
            DimensionMismatch(
                "stoich for reaction $i has length $(length(kr.stoich)), expected $n_sp",
            ),
        )
        S_kin[i, :] = kr.stoich
    end

    # Unique kinetic species (sorted) — Leal et al. 2015 classification:
    # rate-controlling minerals ∪ AS_CRYSTAL species with non-zero stoich
    stoich_crystal = Int[]
    for kr in kinetic_reactions
        for (j, s) in enumerate(kr.stoich)
            if !iszero(s) && aggregate_state(system.species[j]) == AS_CRYSTAL
                push!(stoich_crystal, j)
            end
        end
    end
    idx_kin_unique = sort(unique([idx_kin; stoich_crystal]))

    return KineticsProblem{
        typeof(system),
        typeof(kinetic_reactions),
        typeof(equilibrium_solver),
        typeof(activity_model),
    }(
        system,
        kinetic_reactions,
        initial_state,
        (Float64(safe_ustrip(us"s", tspan[1])), Float64(safe_ustrip(us"s", tspan[2]))),
        activity_model,
        equilibrium_solver,
        S_kin,
        idx_kin,
        idx_kin_unique,
    )
end

"""
    KineticsProblem(system, reactions, initial_state, tspan;
                    activity_model = DiluteSolutionModel(),
                    equilibrium_solver = nothing) -> KineticsProblem

Reaction-centric convenience constructor: build a [`KineticsProblem`](@ref) directly
from a list of [`Reaction`](@ref) objects carrying their kinetics in `properties`.

Each `Reaction` in `reactions` must have a `:rate` entry (a [`KineticFunc`](@ref) or
compatible callable). See [`KineticReaction(cs, rxn)`](@ref) for details.

# Examples

```julia
rxn[:rate]         = parrot_killoh(PK_PARAMS_C3S, "C3S")
rxn[:heat_per_mol] = 114_617.0

kp = KineticsProblem(cs, [rxn_C3S, rxn_C2S, rxn_C3A, rxn_C4AF], state0, (0.0, 7 * 86400.0))
```
"""
function KineticsProblem(
        system::ChemicalSystem,
        reactions::AbstractVector{<:AbstractReaction},
        initial_state::ChemicalState,
        tspan::Tuple;
        activity_model::AbstractActivityModel = DiluteSolutionModel(),
        equilibrium_solver = nothing,
    )
    kinetic_rxns = [KineticReaction(system, rxn) for rxn in reactions]
    return KineticsProblem(
        system, kinetic_rxns, initial_state, tspan;
        activity_model, equilibrium_solver,
    )
end

# ── build_u0 ──────────────────────────────────────────────────────────────────

"""
    build_u0(kp::KineticsProblem) -> Vector{Float64}

Extract the initial mole vector for kinetic (mineral) species from `kp.initial_state`.

The returned vector has length `length(kp.idx_kin_unique)` — one entry per unique
kinetic mineral, in the order defined by `kp.idx_kin_unique`.

# Examples

```julia
u0 = build_u0(kp)
# u0[j] = moles of the j-th unique kinetic mineral at t=0
```
"""
function build_u0(kp::KineticsProblem)
    n_mol = ustrip.(us"mol", kp.initial_state.n)
    return Float64[n_mol[i] for i in kp.idx_kin_unique]
end

# ── build_kinetics_params ─────────────────────────────────────────────────────

"""
    build_kinetics_params(kp::KineticsProblem; ϵ=1e-30) -> NamedTuple

Build the immutable parameter `NamedTuple` `p` passed to the ODE function.

Fields in the returned tuple:
  - `T`: temperature [K] (plain `Float64`).
  - `P`: pressure [Pa].
  - `ϵ`: regularisation floor.
  - `lna_fn`: log-activity closure built from `kp.activity_model`.
  - `kin_rxns`: `kp.kinetic_reactions` (compiled [`KineticFunc`](@ref) objects).
  - `S_kin`: stoichiometric submatrix `[n_reactions × n_species]`.
  - `idx_kin`: mineral indices per reaction.
  - `idx_kin_unique`: unique sorted mineral indices — matches the ODE state `u`.
  - `species_index`: `Dict{String,Int}` mapping PHREEQC formula → species index in system.
    Built once at construction; used to create [`StateView`](@ref)s in the ODE loop
    without allocating a new dict at each step.
  - `n_initial_full`: `Vector{Float64}` of initial moles for **all** species
    (length = n_species). Wrapped in a `StateView` for O(1) named access.
  - `n_full`: mutable `Vector{Float64}` holding the full species mole vector
    (updated in-place at every ODE evaluation).
  - `cp_fns`: vector of Cp°(T) callables [J/(mol·K)] (one per species, `nothing` when
    no Cp° data available). Used by [`SemiAdiabaticCalorimeter`](@ref) to compute
    the variable heat capacity `Cp_cell + Σᵢ nᵢ Cpᵢ°(T)` (Lavergne et al. 2018).
  - `eq_solver`: the equilibrium solver (or `nothing`).
  - `state_ref`: `Ref` to the current [`ChemicalState`](@ref) (updated during integration).
  - `h_fns`: vector of `ΔₐH⁰` callables (or `nothing`) per species.
  - `saved_H`, `saved_t`: mutable buffers for total-enthalpy tracking.
  - `rates_buf`: mutable `Vector{Float64}` of length n_reactions; individual rates at
    each Float64 ODE step. Used by calorimeter `extend_ode!`.
"""
function build_kinetics_params(kp::KineticsProblem; ϵ::Float64 = 1.0e-30)
    state = kp.initial_state
    T = temperature(state)
    P = pressure(state)

    T_K = Float64(ustrip(us"K", T))
    P_Pa = Float64(ustrip(us"Pa", P))

    lna_fn = activity_model(kp.system, kp.activity_model)

    # Pre-built species name → index dict (built once; reused as StateView.index).
    # Two keys per species: PHREEQC formula (e.g. "Ca3SiO5") and species symbol
    # (e.g. "C3S") so that both transition_state (uses formula) and parrot_killoh
    # (uses cement symbol) resolve correctly.
    species_index = Dict{String, Int}()
    for (i, sp) in enumerate(kp.system.species)
        species_index[phreeqc(formula(sp))] = i
        sym = ChemistryLab.symbol(sp)
        if !isempty(sym)
            species_index[sym] = i
        end
    end

    # Full initial mole vector for all species (for n_initial StateView)
    n_initial_full = Float64[ustrip(us"mol", kp.initial_state.n[i]) for i in eachindex(kp.system.species)]

    # Mutable full-system mole vector (updated in-place at each ODE step)
    n_full = Float64[ustrip(us"mol", state.n[i]) for i in eachindex(kp.system.species)]

    # Cp°(T) callables per species [J/(mol·K)], nothing if unavailable
    cp_fns = [haskey(sp, :Cp⁰) ? sp[:Cp⁰] : nothing for sp in kp.system.species]

    # Enthalpy callables for total-enthalpy tracking
    h_fns = [haskey(sp, :ΔₐH⁰) ? sp[:ΔₐH⁰] : nothing for sp in kp.system.species]

    saved_H = Float64[]
    saved_t = Float64[]

    rates_buf = zeros(Float64, length(kp.kinetic_reactions))

    return (
        T = T_K,
        P = P_Pa,
        ϵ = ϵ,
        lna_fn = lna_fn,
        kin_rxns = kp.kinetic_reactions,
        S_kin = kp.S_kin,
        idx_kin = kp.idx_kin,
        idx_kin_unique = kp.idx_kin_unique,
        species_index = species_index,
        n_initial_full = n_initial_full,
        n_full = n_full,
        cp_fns = cp_fns,
        eq_solver = kp.equilibrium_solver,
        state_ref = Ref{ChemicalState}(state),
        h_fns = h_fns,
        saved_H = saved_H,
        saved_t = saved_t,
        rates_buf = rates_buf,
    )
end

# ── build_kinetics_ode ────────────────────────────────────────────────────────

"""
    build_kinetics_ode(kp::KineticsProblem) -> Function

Build the ODE right-hand-side function `f!(du, u, p, t)` for the kinetics problem.

The returned closure implements the Reaktoro-style coupled kinetics–equilibrium
algorithm (Leal et al. 2017):

1. **Update full mole vector**: inject `u` (unique kinetic minerals) into `p.n_full`.
2. **Re-equilibrate** (if `p.eq_solver ≢ nothing`): run aqueous speciation.
3. **Compute log-activities** `lna` via `p.lna_fn`.
4. **Build `StateView`s** (O(1) named access, no dict allocation per step):
   `n_sv`, `lna_sv`, `n_initial_sv` wrapping `p.species_index`.
5. **Evaluate rates**: `r = kr.rate_fn(p.T, p.P, t, n_sv, lna_sv, n_initial_sv)`.
6. **Accumulate `du`**: `du[k] += kr.stoich[sp_idx] * r` for each tracked species.

# Returns

A mutating function `f!(du, u, p, t)` suitable for `ODEProblem`.
"""
function build_kinetics_ode(::KineticsProblem)
    function f!(du, u, p, t)
        T_elt = eltype(u)

        # ── 1. Build full mole vector ──────────────────────────────────────
        if T_elt === Float64
            for (j, idx) in enumerate(p.idx_kin_unique)
                p.n_full[idx] = max(u[j], p.ϵ)
            end
            n_full = p.n_full
        else
            n_full = T_elt.(p.n_full)
            for (j, idx) in enumerate(p.idx_kin_unique)
                n_full[idx] = max(u[j], p.ϵ)
            end
        end

        # ── 2. Re-equilibrate aqueous speciation (Float64 path only) ───────
        if !isnothing(p.eq_solver) && T_elt === Float64
            curr_state = p.state_ref[]
            new_state = ChemicalState(
                curr_state.system,
                p.n_full .* u"mol",
                temperature(curr_state),
                pressure(curr_state),
            )
            try
                eq_result = equilibrate(new_state, p.eq_solver)
                n_eq = ustrip.(us"mol", eq_result.n)
                for i in eachindex(n_eq)
                    if i ∉ p.idx_kin_unique
                        p.n_full[i] = n_eq[i]
                    end
                end
                p.state_ref[] = eq_result
            catch
            end
        end

        # ── 3. Compute log-activities ─────────────────────────────────────
        lna = p.lna_fn(n_full, p)

        # ── 4. Build StateViews (O(1) named access, no dict alloc) ────────
        n_sv = StateView(n_full, p.species_index)
        lna_sv = StateView(lna, p.species_index)
        n_initial_sv = StateView(p.n_initial_full, p.species_index)

        # ── 5–6. Evaluate rates and accumulate ────────────────────────────
        fill!(du, zero(eltype(du)))

        for (i, kr) in enumerate(p.kin_rxns)
            r = kr.rate_fn(p.T, p.P, t, n_sv, lna_sv, n_initial_sv)

            for (k, sp_idx) in enumerate(p.idx_kin_unique)
                s = kr.stoich[sp_idx]
                iszero(s) && continue
                du[k] += s * r
            end

            if T_elt === Float64
                p.rates_buf[i] = r
            end
        end

        return nothing
    end

    return f!
end
