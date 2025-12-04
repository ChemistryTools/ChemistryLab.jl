@kwdef struct EquilibriumProblem{μType<:Function,T<:Number, P, S}
    b::Vector{T}
    A::Matrix{T}
    μ::μType
    u0::Vector{T}
    p::P = SciMLBase.NullParameters()
    lb::Vector{T}
    ub::Vector{T}
    indep_components::Vector{S} = String[]
    dep_components::Vector{S} = String[]
    ϵ::T = 1.e-16
    function EquilibriumProblem(A, μ, u0; p=SciMLBase.NullParameters(),
                                lb=zero(u0), ub=10*sum(u0)*one.(u0),
                                indep_components=fill("", length(b)),
                                dep_components=fill("", length(u0)),
                                ϵ=1.e-16
                                )
        return new{typeof(μ), eltype(u0), typeof(p), eltype(indep_components)}(A*u0, A, μ, u0, p, max.(lb, ϵ), max.(ub, ϵ), indep_components, dep_components, ϵ)
    end
end

function SciMLBase.OptimizationProblem(ep::EquilibriumProblem; kwargs...)
    Gibbs_energy(x, p) = x ⋅ ep.μ(x, p)
    diff_Gibbs_energy!(g, x, p) = g .= ep.μ(x, p)
    cons(res, x, _) = res .= ep.A * x - ep.b

    optf = OptimizationFunction(
        Gibbs_energy,
        Optimization.AutoForwardDiff();
        # grad = diff_Gibbs_energy!,
        cons = cons
    )

    return OptimizationProblem(
        optf,
        ep.u0,
        ep.p,
        lb=ep.lb,
        ub=ep.ub,
        lcons=zeros(size(ep.A, 1)),
        ucons=zeros(size(ep.A, 1)),
        kwargs...
    )
end

SciMLBase.solve(ep::EquilibriumProblem, solver; kwargs...) = solve(OptimizationProblem(ep), solver; kwargs...)
