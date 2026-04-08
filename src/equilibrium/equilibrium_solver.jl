# ── EquilibriumSolver ─────────────────────────────────────────────────────────

"""
    struct EquilibriumSolver{F<:Function, S, V<:Val}

Encapsulates all fixed ingredients of a chemical equilibrium calculation:
the potential function, the SciML solver, and the variable space.

Construct once, call repeatedly with different `ChemicalState` inputs.

# Fields

  - `μ`: chemical potential closure `μ(n, p) -> Vector{Float64}`.
  - `solver`: any Optimization.jl-compatible solver (e.g. `IpoptOptimizer()`).
  - `variable_space`: variable space — `Val(:linear)` or `Val(:log)`.
  - `kwargs`: solver keyword arguments forwarded to `solve`.

# Examples
```julia
julia> cs = ChemicalSystem([
           Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("H+";  aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
           Species("OH-"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
       ]);

julia> solver = EquilibriumSolver(cs, DiluteSolutionModel(), IpoptOptimizer());

julia> solver isa EquilibriumSolver
true
```
"""
struct EquilibriumSolver{F <: Function, S, V <: Val}
    μ::F                    # potential closure — built once from cs and model
    solver::S               # SciML-compatible solver
    variable_space::V       # Val(:linear) or Val(:log)
    kwargs::Base.Pairs      # forwarded to solve
end

"""
    EquilibriumSolver(cs, model, solver; variable_space=Val(:linear), kwargs...)

Construct an `EquilibriumSolver` from a `ChemicalSystem`, an activity model,
and a SciML solver.

The potential function is built once at construction time from `cs` and `model`.
Repeated calls to `solve` with different `ChemicalState` inputs reuse it.

# Arguments

  - `cs`: the `ChemicalSystem` defining species and conservation matrix.
  - `model`: an `AbstractActivityModel` (e.g. `DiluteSolutionModel()`).
  - `solver`: any Optimization.jl solver.
  - `variable_space`: `Val(:linear)` (default) or `Val(:log)`.
  - `kwargs...`: forwarded to the underlying `solve` call (tolerances, verbosity...).
"""
function EquilibriumSolver(
        cs::ChemicalSystem,
        model::AbstractActivityModel,
        solver::S;
        variable_space::V = Val(:linear),
        kwargs...,
    ) where {S, V <: Val}
    μ = build_potentials(cs, model)     # built once — captures indices and constants
    return EquilibriumSolver{typeof(μ), S, V}(μ, solver, variable_space, kwargs)
end

# ── Internal helper: build p from ChemicalState ───────────────────────────────

"""
    _build_params(state::ChemicalState; ϵ=1e-16) -> NamedTuple

Extract dimensionless parameters from a `ChemicalState`.
`ΔₐG⁰overT` is evaluated at the current `T` and `P` of the state.
Units are stripped — compatible with ForwardDiff dual numbers.

# Returned fields

  - `ΔₐG⁰overT`: vector of standard Gibbs energies divided by RT (dimensionless).
  - `T`: temperature in K (plain number, Dual-safe).
  - `P`: pressure in Pa (plain number, Dual-safe).
  - `ϵ`: regularisation floor (default `1e-16`).

`T` and `P` are included so that temperature-dependent activity models
(e.g. [`HKFActivityModel`](@ref) with `temperature_dependent=true`) can
recompute their parameters inside the potential closure.
"""
function _build_params(state::ChemicalState; ϵ::Float64 = 1.0e-16)
    T = temperature(state)
    P = pressure(state)
    R = Constants.R
    RT = R * T                  # keeps units — division below strips them

    # ustrip without forced Float64 conversion — preserves Dual if T is Dual
    ΔₐG⁰overT = [
        ustrip(s[:ΔₐG⁰](T = T, P = P; unit = true) / RT)
            for s in state.system.species
    ]

    T_K = ustrip(us"K", T)   # Quantity{Dual} → Dual, Float64 → Float64
    P_Pa = ustrip(us"Pa", P)

    return (ΔₐG⁰overT = ΔₐG⁰overT, T = T_K, P = P_Pa, ϵ = ϵ)
end

"""
    _build_n0(state::ChemicalState) -> Vector

Extract the dimensionless mole vector from a `ChemicalState`.
Type is inferred from the state — compatible with ForwardDiff dual numbers.
"""
function _build_n0(state::ChemicalState)
    return ustrip.(us"mol", state.n)    # no Float64 cast — type follows eltype(state.n)
end

# ── solve(EquilibriumSolver, ChemicalState) ───────────────────────────────────

"""
    solve(esolver::EquilibriumSolver, state::ChemicalState; ϵ=1e-16) -> ChemicalState

Solve a chemical equilibrium problem from an initial `ChemicalState`.

The conservation matrix `A` is taken from `state.system.SM.A`.
The initial mole vector `n0` and thermodynamic parameters `ΔₐG⁰/RT`
are extracted from `state` at its current `T` and `P`.

Returns a new `ChemicalState` with equilibrium mole amounts,
sharing the same `ChemicalSystem` as the input.

# Arguments

  - `esolver`: the `EquilibriumSolver` to use.
  - `state`: initial state — defines `T`, `P`, and initial composition.
  - `ϵ`: regularization floor (default: `1e-16`).

# Examples
```julia
solver = EquilibriumSolver(cs, DiluteSolutionModel(), IpoptOptimizer();
                           variable_space=Val(:log), abstol=1e-10)
state0 = ChemicalState(cs, n0; T=298.15u"K", P=1u"bar")
state_eq = solve(solver, state0)
```
"""
function SciMLBase.solve(
        esolver::EquilibriumSolver,
        state::ChemicalState;
        ϵ::Float64 = 1.0e-16,
    )
    cs = state.system

    # Extract dimensionless inputs from state
    n0 = _build_n0(state)
    p = _build_params(state; ϵ = ϵ)

    # Ensure initial guess is positive (critical for log variables)
    n0 = max.(n0, ϵ)

    # Build and solve the optimization problem
    prob = EquilibriumProblem(cs.SM.A, esolver.μ, n0; p = p)
    sol = solve(prob, esolver.solver; variable_space = esolver.variable_space, esolver.kwargs...)

    # Wrap result back into a ChemicalState — same system, same T and P
    state_eq = copy(state)                      # shares ChemicalSystem, copies n/T/P
    for (i, nᵢ) in enumerate(sol.u)
        state_eq.n[i] = max(nᵢ, ϵ) * u"mol"   # reattach units, enforce floor
    end
    _update_derived!(state_eq)                  # recompute phases, pH, porosity...

    return state_eq
end

"""
    equilibrate(state::ChemicalState;
                model    = DiluteSolutionModel(),
                solver   = IpoptOptimizer(...),
                variable_space  = Val(:log),
                ϵ        = 1e-16,
                kwargs...) -> ChemicalState

Compute the chemical equilibrium state from an initial `ChemicalState`.

This is a high-level convenience function that wraps `EquilibriumSolver` and
`solve` with sensible defaults. It is intended for users who do not need to
fine-tune the solver.

# Arguments

  - `state`: initial `ChemicalState` — defines the system, T, P, and composition.
  - `model`: activity model (default: `DiluteSolutionModel()`).
  - `solver`: SciML-compatible solver (default: `IpoptOptimizer` with tight tolerances).
  - `variable_space`: variable space — `Val(:linear)` or `Val(:log)` (default).
    `Val(:log)` is more robust for systems spanning many orders of magnitude.
  - `ϵ`: regularization floor for mole amounts (default: `1e-16`).
  - `kwargs...`: additional keyword arguments forwarded to the solver.

# Returns

A new `ChemicalState` at thermodynamic equilibrium, with all derived quantities
(pH, pOH, volumes, porosity, saturation) already computed.

# Examples
```julia
# Minimal usage — all defaults
state_eq = equilibrate(state)

# Custom temperature-dependent run
set_temperature!(state, 350.0u"K")
state_eq = equilibrate(state)

# Custom solver tolerance
state_eq = equilibrate(state; abstol=1e-12, reltol=1e-12)

# Custom activity model (future)
state_eq = equilibrate(state; model=DebyeHuckelModel(A=0.51, B=3.28))
```
"""
function equilibrate(
        state::ChemicalState;
        model::AbstractActivityModel = DiluteSolutionModel(),
        solver = _default_ipopt_solver(),
        variable_space::Val = Val(:linear),
        ϵ::Float64 = 1.0e-16,
        kwargs...,
    )
    esolver = EquilibriumSolver(
        state.system, model, solver;
        variable_space = variable_space,
        kwargs...,
    )
    return solve(esolver, state; ϵ = ϵ)
end

"""
    _default_ipopt_solver() -> IpoptOptimizer

Return an `IpoptOptimizer` instance with tight tolerances suitable for
chemical equilibrium calculations.
"""
function _default_ipopt_solver()
    return IpoptOptimizer(
        acceptable_tol = 1.0e-12,
        dual_inf_tol = 1.0e-12,
        acceptable_iter = 1000,
        constr_viol_tol = 1.0e-12,
        warm_start_init_point = "yes",
    )
end
