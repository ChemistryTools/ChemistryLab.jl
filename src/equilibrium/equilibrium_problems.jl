"""
    EquilibriumProblem

Definition of a chemical equilibrium problem.

# Fields

  - `b`: conservation vector (elemental abundances).
  - `A`: stoichiometric matrix (conservation matrix).
  - `μ`: chemical potential function `μ(n, p)`.
  - `u0`: initial guess for species amounts.
  - `p`: coefficients for the potential function (default: NullParameters).
  - `lb`: lower bounds for species amounts.
  - `ub`: upper bounds for species amounts.

The problem solves for the species distribution that minimizes the Gibbs energy
subject to mass conservation constraints `A * n = b`.
"""
struct EquilibriumProblem{F<:Function, Tb, TA, Tu, P}
    b::Vector{Tb}
    A::Matrix{TA}
    μ::F
    u0::Vector{Tu}
    p::P
    lb::Vector{Tu}
    ub::Vector{Tu}
end

"""
    EquilibriumProblem(A, μ, u0; b=A*u0, p=NullParameters(), lb=fill(Tu(1e-16), length(u0)), ub=maximum(abs.(A))/minimum(abs.(A[.!iszero.(A)]))*sum(u0)*one.(u0))

Construct an `EquilibriumProblem` with the given stoichiometric matrix `A`, chemical potential function `μ`, and initial guess `u0`.

# Arguments

  - `A`: stoichiometric matrix (conservation matrix).
  - `μ`: chemical potential function `μ(n, p)`.
  - `u0`: initial guess for species amounts.
  - `b`: conservation vector (elemental abundances). Defaults to `A * u0`.
  - `p`: coefficients for the potential function. Defaults to `NullParameters()`.
  - `lb`: lower bounds for species amounts. Defaults to `fill(Tu(1e-16), length(u0))`.
  - `ub`: upper bounds for species amounts. Defaults to `maximum(abs.(A))/minimum(abs.(A[.!iszero.(A)]))*sum(u0)*one.(u0)`.

# Returns

An `EquilibriumProblem` instance.
"""
function EquilibriumProblem(
    A::AbstractMatrix{TA},
    μ::F,
    u0::AbstractVector{Tu};
    b::AbstractVector  = A * u0,
    p                  = SciMLBase.NullParameters(),
    lb::AbstractVector = fill(Tu(1e-16), length(u0)),
    ub::AbstractVector = maximum(abs.(A)) / minimum(abs.(A[.!iszero.(A)])) * sum(u0) * one.(u0),
) where {Tu<:Number, TA<:Number, F<:Function}
    ϵ  = 1e-16
    lb = max.(lb, ϵ)
    ub = max.(ub, ϵ)
    # Ensure u0 has no zeros or negative values
    u0 = max.(u0, ϵ)
    Tb = eltype(b)
    return EquilibriumProblem{F, Tb, TA, Tu, typeof(p)}(
        Vector{Tb}(b), Matrix{TA}(A), μ, Vector{Tu}(u0), p,
        Vector{Tu}(lb), Vector{Tu}(ub),
    )
end

"""
    SciMLBase.OptimizationProblem(ep::EquilibriumProblem, variable_space::Val{:linear}; kwargs...)

Convert an `EquilibriumProblem` to an `OptimizationProblem` for linear variables.

# Arguments

  - `ep`: an `EquilibriumProblem` instance.
  - `variable_space`: variable space (`Val(:linear)` for direct mole optimization).
  - `kwargs`: additional keyword arguments to pass to `OptimizationProblem`.

# Returns

An `OptimizationProblem` instance optimized for linear variable space.
"""
function SciMLBase.OptimizationProblem(ep::EquilibriumProblem, ::Val{:linear}; kwargs...)
    # Objective function: minimize Gibbs energy G = n⋅μ(n)
    Gibbs_energy(x, p) = dot(x, ep.μ(x, p))

    # Conservation constraints: A*n = b
    conservation_constraints(res, x, _) = mul!(res, ep.A, x) .-= ep.b

    optf = OptimizationFunction(
        Gibbs_energy,
        Optimization.AutoForwardDiff();
        cons=conservation_constraints,
    )

    return OptimizationProblem(
        optf,
        ep.u0,          # Initial guess in linear space
        ep.p;
        lb=ep.lb,        # Lower bounds (already ensured positive)
        ub=ep.ub,        # Upper bounds (already ensured positive)
        lcons=zeros(size(ep.A, 1)),  # Equality constraints: A*n = b
        ucons=zeros(size(ep.A, 1)),
        kwargs...,
    )
end

"""
    SciMLBase.OptimizationProblem(ep::EquilibriumProblem, variable_space::Val{:log}; kwargs...)

Convert an `EquilibriumProblem` to an `OptimizationProblem` for logarithmic variables.

# Arguments

  - `ep`: an `EquilibriumProblem` instance.
  - `variable_space`: variable space (`Val(:log)` for log(mole) optimization).
  - `kwargs`: additional keyword arguments to pass to `OptimizationProblem`.

# Returns

An `OptimizationProblem` instance optimized for logarithmic variable space.

# Notes

The logarithmic formulation is more robust for systems spanning many orders of magnitude
and automatically enforces positivity of mole amounts through the exp() transformation.
"""
function SciMLBase.OptimizationProblem(ep::EquilibriumProblem, ::Val{:log}; kwargs...)
    # Objective function: minimize Gibbs energy in log-space
    # G = exp(x)⋅μ(exp(x)) where x = log(n)
    Gibbs_energy_log(x, p) = begin
        n = exp.(x)  # Transform to physical space
        return dot(n, ep.μ(n, p))
    end

    # Conservation constraints in log-space: A*exp(x) = b
    conservation_constraints_log(res, x, _) = begin
        n = exp.(x)  # Transform to physical space
        mul!(res, ep.A, n)
        res .-= ep.b
    end

    optf = OptimizationFunction(
        Gibbs_energy_log,
        Optimization.AutoForwardDiff();
        cons=conservation_constraints_log,
    )

    # Transform bounds to log-space
    log_lb = log.(ep.lb)
    log_ub = log.(ep.ub)
    log_u0 = log.(ep.u0)

    return OptimizationProblem(
        optf,
        log_u0,         # Initial guess in log-space
        ep.p;
        lb=log_lb,       # Lower bounds in log-space
        ub=log_ub,       # Upper bounds in log-space
        lcons=zeros(size(ep.A, 1)),  # Equality constraints: A*exp(x) = b
        ucons=zeros(size(ep.A, 1)),
        kwargs...,
    )
end

# Solution transformation methods using multiple dispatch
_solution_transform(::Val{:linear}) = identity
_solution_transform(::Val{:log}) = exp

"""
    SciMLBase.solve(ep::EquilibriumProblem, solver; variable_space=Val(:linear), kwargs...)

Solve an `EquilibriumProblem` using the specified solver and variable space.

# Arguments

  - `ep`: an `EquilibriumProblem` instance.
  - `solver`: the solver to use.
  - `variable_space`: the variable space to use (`Val(:linear)` or `Val(:log)`). Defaults to `Val(:linear)`.
  - `kwargs`: additional keyword arguments to pass to `solve`.

# Returns

The solution to the `EquilibriumProblem` with variables transformed back to physical space.
"""
function SciMLBase.solve(ep::EquilibriumProblem, solver; variable_space=Val(:linear), kwargs...)
    # Create and solve the optimization problem
    opt_prob = OptimizationProblem(ep, variable_space)
    sol = solve(opt_prob, solver; kwargs...)

    # Transform solution back to physical space using multiple dispatch
    transform_fn = _solution_transform(variable_space)
    sol.u .= transform_fn.(sol.u)
    return sol
end
