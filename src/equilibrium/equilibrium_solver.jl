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

# ── equilibrate(ChemicalState) ────────────────────────────────────────────────

"""
    equilibrate(state::ChemicalState;
                model          = DiluteSolutionModel(),
                solver         = IpoptOptimizer(...),
                variable_space = Val(:linear),
                ϵ              = 1e-16,
                kwargs...) -> ChemicalState

Compute the chemical equilibrium state from an initial `ChemicalState`.

Requires `Optimization` and `OptimizationIpopt` to be loaded:

```julia
using Optimization, OptimizationIpopt
using ChemistryLab
state_eq = equilibrate(state)
```

# Arguments

  - `state`: initial `ChemicalState` — defines the system, T, P, and composition.
  - `model`: activity model (default: `DiluteSolutionModel()`).
  - `solver`: SciML-compatible solver (default: `IpoptOptimizer` with tight tolerances).
  - `variable_space`: `Val(:linear)` (default) or `Val(:log)`.
  - `ϵ`: regularization floor for mole amounts (default: `1e-16`).
  - `kwargs...`: additional keyword arguments forwarded to the solver.
"""
function equilibrate(args...; kwargs...)
    return error(
        "equilibrate requires Optimization.jl and OptimizationIpopt.jl to be loaded.\n" *
        "Add `using Optimization, OptimizationIpopt` before calling this function.",
    )
end

# ── solve(EquilibriumSolver{OptimaOptimizer}, ChemicalState) ──────────────────
# Available whenever OptimaSolver is loaded (direct dependency).
# Uses SciMLBase.OptimizationFunction with NoAD — OptimaOptimizer handles
# gradients internally via ForwardDiff.

function _build_optima_opt_prob(ep::EquilibriumProblem, μ, ::Val{:linear})
    f_gibbs(x, q) = LinearAlgebra.dot(x, μ(x, q))
    cons!(res, x, _) = LinearAlgebra.mul!(res, ep.A, x) .-= ep.b
    optf = SciMLBase.OptimizationFunction{true}(f_gibbs; cons = cons!)
    return SciMLBase.OptimizationProblem(
        optf, ep.u0, ep.p;
        lb = ep.lb, ub = ep.ub,
        lcons = zeros(size(ep.A, 1)),
        ucons = zeros(size(ep.A, 1)),
    )
end

function _build_optima_opt_prob(ep::EquilibriumProblem, μ, ::Val{:log})
    f_gibbs(x, q) = (n = exp.(x); LinearAlgebra.dot(n, μ(n, q)))
    cons!(res, x, _) = (n = exp.(x); LinearAlgebra.mul!(res, ep.A, n); res .-= ep.b)
    optf = SciMLBase.OptimizationFunction{true}(f_gibbs; cons = cons!)
    return SciMLBase.OptimizationProblem(
        optf, log.(ep.u0), ep.p;
        lb = log.(ep.lb), ub = log.(ep.ub),
        lcons = zeros(size(ep.A, 1)),
        ucons = zeros(size(ep.A, 1)),
    )
end

"""
    SciMLBase.solve(esolver::EquilibriumSolver{F,<:OptimaOptimizer,V},
                   state::ChemicalState; ϵ=1e-16) -> ChemicalState

Solve a chemical equilibrium problem using an `OptimaOptimizer` solver.
Available without loading `Optimization.jl` / `OptimizationIpopt.jl`.
"""
function SciMLBase.solve(
        esolver::EquilibriumSolver{F, <:OptimaOptimizer, V},
        state::ChemicalState;
        ϵ::Float64 = 1.0e-16,
    ) where {F, V}
    n0 = max.(_build_n0(state), ϵ)
    p = _build_params(state; ϵ = ϵ)

    prob = EquilibriumProblem(state.system.SM.A, esolver.μ, n0; p = p)
    opt_prob = _build_optima_opt_prob(prob, esolver.μ, esolver.variable_space)

    sol = SciMLBase.solve(opt_prob, esolver.solver; esolver.kwargs...)
    transform = _solution_transform(esolver.variable_space)

    state_eq = copy(state)
    for (i, nᵢ) in enumerate(sol.u)
        state_eq.n[i] = max(transform(nᵢ), ϵ) * u"mol"
    end
    _update_derived!(state_eq)

    return state_eq
end
