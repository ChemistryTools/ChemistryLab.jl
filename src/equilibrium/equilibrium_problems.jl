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
  - `indep_components`: names/symbols of independent components.
  - `dep_components`: names/symbols of dependent components (species).

The problem solves for the species distribution that minimizes the Gibbs energy
subject to mass conservation constraints `A * n = b`.
"""
@kwdef struct EquilibriumProblem{μType<:Function,T<:Number,P,S}
    b::Vector{T}
    A::Matrix{T}
    μ::μType
    u0::Vector{T}
    p::P = SciMLBase.NullParameters()
    lb::Vector{T}
    ub::Vector{T}
    indep_components::Vector{S} = String[]
    dep_components::Vector{S} = String[]
    function EquilibriumProblem(
        A,
        μ,
        u0;
        p=SciMLBase.NullParameters(),
        lb=zero(u0),
        ub=maximum(abs.(A)) / minimum(abs.(A[.!iszero.(A)])) * sum(u0) * one.(u0),
        indep_components=fill("", length(b)),
        dep_components=fill("", length(u0)),
    )
        ϵ = 1.e-16
        return new{typeof(μ),eltype(u0),typeof(p),eltype(indep_components)}(
            A * u0, A, μ, u0, p, max.(lb, ϵ), max.(ub, ϵ), indep_components, dep_components
        )
    end
end

function SciMLBase.OptimizationProblem(ep::EquilibriumProblem, ::Val{:linear}; kwargs...)
    Gibbs_energy(x, p) = x ⋅ ep.μ(x, p)
    diff_Gibbs_energy!(g, x, p) = g .= ep.μ(x, p)
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

function SciMLBase.solve(ep::EquilibriumProblem, vartype, solver; kwargs...)
    return solve(OptimizationProblem(ep, vartype), solver; kwargs...)
end
