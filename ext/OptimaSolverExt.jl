module OptimaSolverExt

using ChemistryLab
import ChemistryLab:
    EquilibriumProblem,
    EquilibriumSolver,
    ChemicalState,
    _build_params,
    _build_n0,
    _solution_transform,
    _update_derived!
using OptimaSolver: OptimaOptimizer
using SciMLBase
using LinearAlgebra: dot, mul!
using DynamicQuantities

# ── OptimizationProblem helpers (NoAD — OptimaOptimizer handles gradients) ────

function _build_optima_opt_prob(ep::EquilibriumProblem, μ, ::Val{:linear})
    f_gibbs(x, q) = dot(x, μ(x, q))
    cons!(res, x, _) = mul!(res, ep.A, x) .-= ep.b
    optf = SciMLBase.OptimizationFunction{true}(f_gibbs; cons = cons!)
    return SciMLBase.OptimizationProblem(
        optf, ep.u0, ep.p;
        lb = ep.lb, ub = ep.ub,
        lcons = zeros(size(ep.A, 1)),
        ucons = zeros(size(ep.A, 1)),
    )
end

function _build_optima_opt_prob(ep::EquilibriumProblem, μ, ::Val{:log})
    f_gibbs(x, q) = (n = exp.(x); dot(n, μ(n, q)))
    cons!(res, x, _) = (n = exp.(x); mul!(res, ep.A, n); res .-= ep.b)
    optf = SciMLBase.OptimizationFunction{true}(f_gibbs; cons = cons!)
    return SciMLBase.OptimizationProblem(
        optf, log.(ep.u0), ep.p;
        lb = log.(ep.lb), ub = log.(ep.ub),
        lcons = zeros(size(ep.A, 1)),
        ucons = zeros(size(ep.A, 1)),
    )
end

# ── solve(EquilibriumSolver{OptimaOptimizer}, ChemicalState) ──────────────────

"""
    SciMLBase.solve(esolver::EquilibriumSolver{F,<:OptimaOptimizer,V},
                   state::ChemicalState; ϵ=1e-16) -> ChemicalState

Solve a chemical equilibrium problem using an `OptimaOptimizer` solver.
Loaded automatically when `using OptimaSolver` is active.
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

end # module OptimaSolverExt
