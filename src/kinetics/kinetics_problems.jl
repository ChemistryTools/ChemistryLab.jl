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
  - `kinetic_reactions`: vector of [`KineticReaction`](@ref) objects (one per
    kinetically controlled mineral).
  - `initial_state`: [`ChemicalState`](@ref) providing initial mole amounts, T, P.
  - `tspan`: `(t_start, t_end)` time interval [s].
  - `activity_model`: [`AbstractActivityModel`](@ref) used at every ODE step to
    compute log-activities (and thus saturation ratios).
  - `equilibrium_solver`: an [`EquilibriumSolver`](@ref) re-equilibrating the fast
    (aqueous) species at every ODE evaluation. `nothing` disables re-speciation.
  - `S_kin`: stoichiometric submatrix `[n_kin × n_species]` — row `i` is the
    stoichiometric coefficient vector for kinetic reaction `i`.
  - `idx_kin`: integer indices of the kinetic mineral species in `system`.

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
    equilibrium_solver::ES   # EquilibriumSolver or nothing
    S_kin::Matrix{Float64}   # [n_kin × n_species]
    idx_kin::Vector{Int}     # indices of kinetic minerals in system
end

"""
    KineticsProblem(system, kinetic_reactions, initial_state, tspan;
                    activity_model = DiluteSolutionModel(),
                    equilibrium_solver = nothing) -> KineticsProblem

Construct a [`KineticsProblem`](@ref).

The stoichiometric submatrix `S_kin` and kinetic species index vector `idx_kin`
are computed automatically from `kinetic_reactions`.

# Arguments

  - `system`: a [`ChemicalSystem`](@ref).
  - `kinetic_reactions`: `AbstractVector` of [`KineticReaction`](@ref).
  - `initial_state`: a [`ChemicalState`](@ref) with initial moles, T, P.
  - `tspan`: `(t0, tf)` time interval in seconds.
  - `activity_model`: activity model for log-activity computation at each step
    (default: [`DiluteSolutionModel`](@ref)).
  - `equilibrium_solver`: an [`EquilibriumSolver`](@ref) for re-speciation of
    the aqueous phase at each ODE function evaluation (Reaktoro-style coupling).
    Pass `nothing` (default) to skip re-speciation.

# Examples

```julia
kp = KineticsProblem(
    cs, [kr_calcite], state0, (0.0, 7200.0);
    activity_model   = HKFActivityModel(),
    equilibrium_solver = EquilibriumSolver(cs, HKFActivityModel(), my_solver),
)
```
"""
function KineticsProblem(
        system::ChemicalSystem,
        kinetic_reactions::AbstractVector,
        initial_state::ChemicalState,
        tspan::Tuple{<:Real, <:Real};
        activity_model::AbstractActivityModel = DiluteSolutionModel(),
        equilibrium_solver = nothing,
    )
    n_sp = length(system.species)

    # Build stoichiometric submatrix and mineral index vector
    idx_kin = Int[kr.idx_mineral for kr in kinetic_reactions]
    n_kin = length(kinetic_reactions)
    S_kin = zeros(Float64, n_kin, n_sp)
    for (i, kr) in enumerate(kinetic_reactions)
        length(kr.stoich) == n_sp || throw(
            DimensionMismatch(
                "stoich for reaction $i has length $(length(kr.stoich)), expected $n_sp"
            )
        )
        S_kin[i, :] = kr.stoich
    end

    return KineticsProblem{
        typeof(system),
        typeof(kinetic_reactions),
        typeof(equilibrium_solver),
        typeof(activity_model),
    }(
        system,
        kinetic_reactions,
        initial_state,
        (Float64(tspan[1]), Float64(tspan[2])),
        activity_model,
        equilibrium_solver,
        S_kin,
        idx_kin,
    )
end

# ── build_u0 ──────────────────────────────────────────────────────────────────

"""
    build_u0(kp::KineticsProblem) -> Vector{Float64}

Extract the initial mole vector for kinetic (mineral) species from `kp.initial_state`.

The returned vector has length `length(kp.idx_kin)`, one entry per kinetic mineral,
in the order defined by `kp.kinetic_reactions`.

# Examples

```julia
u0 = build_u0(kp)
# u0[i] = moles of kinetic mineral i at t=0
```
"""
function build_u0(kp::KineticsProblem)
    n_mol = ustrip.(us"mol", kp.initial_state.n)
    return Float64[n_mol[i] for i in kp.idx_kin]
end

# ── build_kinetics_params ─────────────────────────────────────────────────────

"""
    build_kinetics_params(kp::KineticsProblem; ϵ=1e-30) -> NamedTuple

Build the immutable parameter `NamedTuple` `p` passed to the ODE function.

Fields in the returned tuple (accessible as `p.T`, `p.ΔₐG⁰overT`, etc.):
  - `T`: temperature [K] (plain `Float64`, Dual-safe).
  - `P`: pressure [Pa].
  - `ΔₐG⁰overT`: standard Gibbs energies divided by RT (length = n_species).
  - `ϵ`: regularisation floor.
  - `lna_fn`: log-activity closure built from `kp.activity_model`.
  - `kin_rxns`: `kp.kinetic_reactions` (rate models + surface models).
  - `S_kin`: stoichiometric submatrix `[n_kin × n_species]`.
  - `idx_kin`: indices of kinetic minerals.
  - `n_full`: mutable `Vector{Float64}` holding the full species mole vector
    (updated in-place at every ODE evaluation from the current `u` + aqueous state).
  - `eq_solver`: the equilibrium solver (or `nothing`).
  - `state_ref`: a `Ref` to the current best-known [`ChemicalState`](@ref)
    (updated in-place during integration when `eq_solver` is active).
  - `h_fns`: vector of standard-enthalpy callables `sp[:ΔₐH⁰]` (or `nothing` for
    species without enthalpy data). Used by the total-enthalpy `DiscreteCallback` in
    `KineticsOrdinaryDiffEqExt` to capture heat from both kinetic and equilibrium
    reactions.
  - `saved_H`: mutable `Vector{Float64}` filled by the `DiscreteCallback` with the
    total system enthalpy H(t) at each accepted ODE step.
  - `saved_t`: mutable `Vector{Float64}` with the corresponding times.
  - `n_initial`: `Vector{Float64}` of initial moles for each kinetic mineral (length = n_kin).
    Passed to rate models that require it (e.g. [`ParrotKillohRateModel`](@ref)).
"""
function build_kinetics_params(kp::KineticsProblem; ϵ::Float64 = 1.0e-30)
    state = kp.initial_state
    T = temperature(state)
    P = pressure(state)
    R = Constants.R
    RT = R * T

    ΔₐG⁰overT = Float64[
        ustrip(s[:ΔₐG⁰](T = T, P = P; unit = true) / RT)
            for s in state.system.species
    ]

    T_K = ustrip(us"K", T)
    P_Pa = ustrip(us"Pa", P)

    lna_fn = activity_model(kp.system, kp.activity_model)

    # Mutable full-system mole vector (updated in-place at each ODE step)
    n_full = ustrip.(us"mol", state.n)

    # Enthalpy functions for total-enthalpy tracking (supports equilibrium heat)
    # One entry per species: the ΔₐH⁰ callable, or nothing if unavailable.
    h_fns = [haskey(sp, :ΔₐH⁰) ? sp[:ΔₐH⁰] : nothing for sp in state.system.species]

    # Mutable buffers filled by the DiscreteCallback in KineticsOrdinaryDiffEqExt
    saved_H = Float64[]
    saved_t = Float64[]

    # Initial moles per kinetic mineral — needed by models such as ParrotKillohRateModel
    n_mol_full = ustrip.(us"mol", state.n)
    n_initial = Float64[n_mol_full[i] for i in kp.idx_kin]

    return (
        T = T_K,
        P = P_Pa,
        ΔₐG⁰overT = ΔₐG⁰overT,
        ϵ = ϵ,
        lna_fn = lna_fn,
        kin_rxns = kp.kinetic_reactions,
        S_kin = kp.S_kin,
        idx_kin = kp.idx_kin,
        n_full = n_full,
        eq_solver = kp.equilibrium_solver,
        state_ref = Ref{ChemicalState}(state),
        h_fns = h_fns,
        saved_H = saved_H,
        saved_t = saved_t,
        n_initial = n_initial,
    )
end

# ── build_kinetics_ode ────────────────────────────────────────────────────────

"""
    build_kinetics_ode(kp::KineticsProblem) -> Function

Build the ODE right-hand-side function `f!(du, u, p, t)` for the kinetics problem.

The returned closure implements the Reaktoro-style coupled kinetics–equilibrium
algorithm (Leal et al. 2017):

1. **Update full mole vector**: inject the current ODE state `u` (kinetic minerals)
   into the full species mole vector `p.n_full`.
2. **Re-equilibrate** (if `p.eq_solver ≢ nothing`): run
   `equilibrate(state, eq_solver)` on the aqueous phase, updating `p.n_full`
   with the new aqueous speciation. This couples fast equilibrium reactions to
   the slow kinetic ODE.
3. **Compute log-activities** `lna` from `p.n_full` via `p.lna_fn`.
4. **Build catalyst lookup** `lna_dict`: a `Dict{String, Float64}` mapping PHREEQC
   formula strings to their log-activities (for catalyst terms in rate models).
5. **Compute saturation ratio** Ω for each kinetic reaction.
6. **Evaluate rate** r [mol/s] for each kinetic reaction.
7. **Set `du`** = `S_kin[i, :]` inner-product with `r[i]` for the mineral rows
   (i.e. `du[i] = r[i]` since S_kin[i, idx_kin[i]] = -1 for a dissolution reaction).

The state vector `u` has length `n_kin` (number of kinetic mineral species).

**AD note**: the closure uses `eltype(u)` for type promotion, so ForwardDiff
Dual numbers flow through provided `p.eq_solver = nothing` (re-speciation
currently breaks Dual propagation since it calls a nonlinear solver).

# Returns

A mutating function `f!(du, u, p, t)` suitable for `ODEProblem`.

See also: [`KineticsProblem`](@ref), [`build_u0`](@ref), [`build_kinetics_params`](@ref).
"""
function build_kinetics_ode(kp::KineticsProblem)

    # Pre-build species formula → index dict for catalyst lookup
    formula_to_idx = Dict{String, Int}(
        phreeqc(formula(sp)) => i for (i, sp) in enumerate(kp.system.species)
    )

    function f!(du, u, p, ::Any)
        T_elt = eltype(u)

        # ── 1. Build AD-compatible full mole vector ────────────────────────
        # • Float64 path (normal evaluation): update p.n_full in-place — no alloc.
        # • Dual path (Jacobian via ForwardDiff): allocate a typed local copy;
        #   p.n_full cannot hold Dual values (it is Vector{Float64}).
        if T_elt === Float64
            for (j, idx) in enumerate(p.idx_kin)
                p.n_full[idx] = max(u[j], p.ϵ)
            end
            n_full = p.n_full
        else
            n_full = T_elt.(p.n_full)
            for (j, idx) in enumerate(p.idx_kin)
                n_full[idx] = max(u[j], p.ϵ)
            end
        end

        # ── 2. Re-equilibrate aqueous speciation (Reaktoro-style coupling) ──
        # Skipped during Jacobian evaluation (nonlinear solver breaks Dual flow).
        if !isnothing(p.eq_solver) && T_elt === Float64
            # Build a ChemicalState from current n_full, T, P
            curr_state = p.state_ref[]
            new_state = ChemicalState(
                curr_state.system,
                p.n_full .* u"mol",
                temperature(curr_state),
                pressure(curr_state),
            )
            try
                eq_result = equilibrate(new_state, p.eq_solver)
                # Update aqueous portion of n_full only
                n_eq = ustrip.(us"mol", eq_result.n)
                for i in eachindex(n_eq)
                    if i ∉ p.idx_kin
                        p.n_full[i] = n_eq[i]
                    end
                end
                p.state_ref[] = eq_result
            catch
                # If equilibration fails, keep current n_full (best effort)
            end
        end

        # ── 3. Compute log-activities ─────────────────────────────────────
        lna = p.lna_fn(n_full, p)

        # ── 4. Build catalyst log-activity lookup ─────────────────────────
        lna_dict = Dict{String, eltype(lna)}(
            k => lna[v] for (k, v) in formula_to_idx
        )

        # ── 5–6. Compute saturation ratios and rates ──────────────────────
        for (i, kr) in enumerate(p.kin_rxns)
            Ω = saturation_ratio(
                kr.stoich, lna, p.ΔₐG⁰overT; ϵ = p.ϵ
            )
            n_mineral = max(u[i], p.ϵ)
            M_mineral = molar_mass(kr)
            A = surface_area(kr.surface_model, n_mineral, M_mineral)
            r = kr.rate_model(;
                T = p.T,
                Ω = Ω,
                A_surface = A,
                lna_dict = lna_dict,
                ϵ = p.ϵ,
                n_current = n_mineral,
                n_initial = p.n_initial[i],
            )
            # ── 7. du[i] = net rate of change of mineral i [mol/s] ────────
            # By sign convention: r > 0 = dissolution → mineral decreases
            du[i] = -r
        end

        return nothing
    end

    return f!
end
