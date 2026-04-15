using DynamicQuantities
using ForwardDiff
using OrderedCollections
using SciMLBase

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

# ── Default solver factory ────────────────────────────────────────────────────

"""
    _DEFAULT_SOLVER_FACTORY

Internal `Ref{Union{Nothing, Function}}` — populated by extension `__init__`
functions to register a default solver factory.

  - `OptimizationIpoptExt.__init__`: registers only if nothing is set (low priority).
  - `OptimaSolverExt.__init__`: always overrides (high priority).

Result: OptimaSolver wins whenever loaded, regardless of load order.
"""
const _DEFAULT_SOLVER_FACTORY = Ref{Union{Nothing, Function}}(nothing)

# ── equilibrate ──────────────────────────────────────────────────────────────

"""
    equilibrate(state::ChemicalState, solver; model=..., variable_space=..., ϵ=...) -> ChemicalState
    equilibrate(state::ChemicalState; kwargs...) -> ChemicalState

Compute the chemical equilibrium state by minimising the Gibbs free energy.

**Two-argument form** (solver explicit, always available once an extension is loaded):

```julia
using Optimization, OptimizationIpopt
state_eq = equilibrate(state, IpoptOptimizer())

using OptimaSolver
state_eq = equilibrate(state, OptimaOptimizer())
```

**One-argument form** (uses the default solver of the loaded extension):

```julia
state_eq = equilibrate(state)
```

When both extensions are loaded, `OptimaSolverExt` always takes priority.

# Arguments

  - `state`: initial `ChemicalState` — defines the system, T, P, and composition.
  - `solver`: any SciML-compatible solver (e.g. `IpoptOptimizer()`, `OptimaOptimizer()`).
  - `model`: activity model (default: `DiluteSolutionModel()`).
  - `variable_space`: `Val(:linear)` (default) or `Val(:log)`.
  - `ϵ`: regularisation floor for mole amounts (default: `1e-16`).
  - `kwargs...`: forwarded to the underlying solver.
"""
function equilibrate(
        state::ChemicalState,
        solver;
        model::AbstractActivityModel = DiluteSolutionModel(),
        variable_space::Val = Val(:linear),
        ϵ::Float64 = 1.0e-16,
        kwargs...,
    )
    esolver = EquilibriumSolver(
        state.system, model, solver;
        variable_space = variable_space,
        kwargs...,
    )
    return SciMLBase.solve(esolver, state; ϵ = ϵ)
end

function equilibrate(state::ChemicalState; kwargs...)
    f = _DEFAULT_SOLVER_FACTORY[]
    if isnothing(f)
        error(
            "equilibrate without explicit solver requires an extension.\n" *
                "Add `using Optimization, OptimizationIpopt` or `using OptimaSolver`, " *
                "or call `equilibrate(state, solver; ...)` explicitly.",
        )
    end
    return equilibrate(state, f(); kwargs...)
end
