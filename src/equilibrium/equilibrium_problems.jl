using DynamicQuantities

"""
    EquilibriumProblem

Definition of a chemical equilibrium problem.

# Fields

  - `b`: conservation vector (elemental abundances).
  - `A`: stoichiometric matrix (conservation matrix).
  - `μ`: chemical potential function `μ(n, p)`.
  - `u0`: initial guess for species amounts.
  - `p`: coefficients for the potential function (default: `nothing`).
  - `lb`: lower bounds for species amounts.
  - `ub`: upper bounds for species amounts.

The problem solves for the species distribution that minimizes the Gibbs energy
subject to mass conservation constraints `A * n = b`.
"""
struct EquilibriumProblem{F <: Function, Tb, TA, Tu, P}
    b::Vector{Tb}
    A::Matrix{TA}
    μ::F
    u0::Vector{Tu}
    p::P
    lb::Vector{Tu}
    ub::Vector{Tu}
end

"""
    EquilibriumProblem(A, μ, u0; b=A*u0, p=nothing, lb=fill(Tu(1e-16), length(u0)), ub=maximum(abs.(A))/minimum(abs.(A[.!iszero.(A)]))*sum(u0)*one.(u0))

Construct an `EquilibriumProblem` with the given stoichiometric matrix `A`, chemical potential function `μ`, and initial guess `u0`.

# Arguments

  - `A`: stoichiometric matrix (conservation matrix).
  - `μ`: chemical potential function `μ(n, p)`.
  - `u0`: initial guess for species amounts.
  - `b`: conservation vector (elemental abundances). Defaults to `A * u0`.
  - `p`: coefficients for the potential function. Defaults to `nothing`.
  - `lb`: lower bounds for species amounts. Defaults to `fill(Tu(1e-16), length(u0))`.
  - `ub`: upper bounds for species amounts. Defaults to `maximum(abs.(A))/minimum(abs.(A[.!iszero.(A)]))*sum(u0)*one.(u0)`.

# Returns

An `EquilibriumProblem` instance.
"""
function EquilibriumProblem(
        A::AbstractMatrix{TA},
        μ::F,
        u0::AbstractVector{Tu};
        b::AbstractVector = A * u0,
        p = nothing,
        lb::AbstractVector = fill(Tu(1.0e-16), length(u0)),
        ub::AbstractVector = maximum(abs.(A)) / minimum(abs.(A[.!iszero.(A)])) * sum(u0) * one.(u0),
    ) where {Tu <: Number, TA <: Number, F <: Function}
    ϵ = 1.0e-16
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

# Solution transformation methods using multiple dispatch
_solution_transform(::Val{:linear}) = identity
_solution_transform(::Val{:log}) = exp
