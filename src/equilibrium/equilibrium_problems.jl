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
    Tb = eltype(b)
    return EquilibriumProblem{F, Tb, TA, Tu, typeof(p)}(
        Vector{Tb}(b), Matrix{TA}(A), μ, Vector{Tu}(u0), p,
        Vector{Tu}(lb), Vector{Tu}(ub),
    )
end

"""
    SciMLBase.OptimizationProblem(ep::EquilibriumProblem, ::Val{:linear}; kwargs...)

Convert an `EquilibriumProblem` to an `OptimizationProblem` for linear variables.

# Arguments

  - `ep`: an `EquilibriumProblem` instance.
  - `kwargs`: additional keyword arguments to pass to `OptimizationProblem`.

# Returns

An `OptimizationProblem` instance.
"""
function SciMLBase.OptimizationProblem(ep::EquilibriumProblem, ::Val{:linear}; kwargs...)
    Gibbs_energy(x, p) = x ⋅ ep.μ(x, p)
    # diff_Gibbs_energy!(g, x, p) = g .= ep.μ(x, p)
    cons(res, x, _) = res .= ep.A * x - ep.b

    optf = OptimizationFunction(
        Gibbs_energy,
        Optimization.AutoForwardDiff();
        # grad = diff_Gibbs_energy!,
        cons=cons,
    )

    return OptimizationProblem(
        optf,
        ep.u0,
        ep.p;
        lb=ep.lb,
        ub=ep.ub,
        lcons=zeros(size(ep.A, 1)),
        ucons=zeros(size(ep.A, 1)),
        kwargs...,
    )
end

"""
    SciMLBase.OptimizationProblem(ep::EquilibriumProblem, ::Val{:log}; kwargs...)

Convert an `EquilibriumProblem` to an `OptimizationProblem` for logarithmic variables.

# Arguments

  - `ep`: an `EquilibriumProblem` instance.
  - `kwargs`: additional keyword arguments to pass to `OptimizationProblem`.

# Returns

An `OptimizationProblem` instance.
"""
function SciMLBase.OptimizationProblem(ep::EquilibriumProblem, ::Val{:log}; kwargs...)
    Gibbs_energy(x, p) = (y = exp.(x); y ⋅ ep.μ(y, p))
    cons(res, x, _) = res .= ep.A * exp.(x) - ep.b

    optf = OptimizationFunction(Gibbs_energy, Optimization.AutoForwardDiff(); cons=cons)

    return OptimizationProblem(
        optf,
        log.(ep.u0),
        ep.p;
        lb=log.(ep.lb),
        ub=log.(ep.ub),
        lcons=zeros(size(ep.A, 1)),
        ucons=zeros(size(ep.A, 1)),
        kwargs...,
    )
end

"""
    SciMLBase.solve(ep::EquilibriumProblem, solver, vartype=Val(:linear); kwargs...)

Solve an `EquilibriumProblem` using the specified solver and variable type.

# Arguments

  - `ep`: an `EquilibriumProblem` instance.
  - `solver`: the solver to use.
  - `vartype`: the type of variables to use. Defaults to `Val(:linear)`.
  - `kwargs`: additional keyword arguments to pass to `solve`.

# Returns

The solution to the `EquilibriumProblem`.
"""
function SciMLBase.solve(ep::EquilibriumProblem, solver, vartype=Val(:linear); kwargs...)
    return solve(OptimizationProblem(ep, vartype), solver; kwargs...)
end
