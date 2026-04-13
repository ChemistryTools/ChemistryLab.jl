using Pkg
Pkg.activate(@__DIR__)
# Pkg.add(path="../ChemistryLab.jl")

using Revise, ChemistryLab
using Optimization, OptimizationMOI, OptimizationOptimJL, OptimizationIpopt
using Plots
using SparseArrays

R, T = 8.314, 273.15 + 25
RT = R * T

function μ(n, p)
    pp = NamedTuple(p)
    ntot = sum(n)
    nᵋ = max.(n, pp.ϵ)
    δ = log(55.5)
    return pp.ΔₐG⁰ / RT + [i == pp.idxsolvent ? 0 : δ for i in eachindex(n)] + log.(nᵋ ./ ntot)
end

function _n0(p)
    pp = NamedTuple(p)
    ϵ = pp.ϵ
    ca, Va, cb, Vb = pp.ca, pp.Va, pp.cb, pp.Vb
    V = Va + Vb

    n_H2O = 55.5 * V
    n_H = 1.0e-7V
    n_OH = 1.0e-7V

    n_AH = max(ca * Va, ϵ)
    n_A = ϵ
    n_BOH = max(cb * Vb, ϵ)
    n_B = ϵ

    return [n_H2O, n_H, n_OH, n_AH, n_A, n_BOH, n_B]
end

A = [
    2 1 1 1 0 1 0
    1 0 1 0 0 1 0
    0 0 0 1 1 0 0
    0 0 0 0 0 1 1
]

lsp = [:H₂O, :H⁺, :OH⁻, :AH, :A⁻, :BOH, :B⁺]
ΔₐG⁰ = [-237.18e3, 0.0e3, -157.27e3, -396.46e3, -369.31e3, -417.98e3, -261.9e3]
for (sp, g) in zip(lsp, ΔₐG⁰)
    symsp = String(sp)
    @eval begin
        $sp = Species($symsp; aggregate_state = AS_AQUEOUS, class = $symsp == "H₂O" ? SC_AQSOLVENT : SC_AQSOLUTE)
        $sp.ΔₐG⁰ = $g
    end
end
species = eval.(lsp)
SM = StoichMatrix(species)
A, indep_comp, dep_comp = SM.A, SM.primaries, SM.species
indices_in = findall(x -> x in indep_comp, dep_comp)
indices_out = setdiff(eachindex(dep_comp), indices_in)
indices = [indices_out; indices_in]
dep_comp = dep_comp[indices]
A = A[:, indices]
pprint(A, symbol.(indep_comp), symbol.(dep_comp))

A = sparse(A)

function _ub(p)
    pp = NamedTuple(p)
    ϵ = pp.ϵ
    ca, Va, cb, Vb = pp.ca, pp.Va, pp.cb, pp.Vb
    V = Va + Vb

    return [2 * 55.5V, 1.0e-1V, 1.0e-1V, max(ca * Va, ϵ), max(ca * Va, ϵ), max(cb * Vb, ϵ), max(cb * Vb, ϵ)][indices]
end

# Data acide acétique CH₃COOH et soude NaOH
Ka = exp(-(A⁻.ΔₐG⁰ + H⁺.ΔₐG⁰ - AH.ΔₐG⁰) / RT)
pKa = -log10(Ka)
Kb = exp(-(B⁺.ΔₐG⁰ + OH⁻.ΔₐG⁰ - BOH.ΔₐG⁰) / RT)
pKb = -log10(Kb)


p = [
    :ΔₐG⁰ => getproperty.(dep_comp, :ΔₐG⁰),
    :idxsolvent => findfirst(x -> x.class == SC_AQSOLVENT, dep_comp),
    :ϵ => 1.0e-16, :ca => 0.1, :Va => 0.1, :cb => 0.1, :Vb => 0.0,
]

prob = EquilibriumProblem(A, μ, _n0(p)[indices], ub = _ub(p), p = p)
opt = IpoptOptimizer(
    acceptable_tol = 1.0e-8,
    dual_inf_tol = 1.0e-8,
    acceptable_iter = 100,
    constr_viol_tol = 1.0e-8,
    # compl_inf_tol = 1e-4,
    # mu_strategy = "adaptive",
    warm_start_init_point = "yes",
    # expect_infeasible_problem = "yes",
)
sol = exp.(solve(prob, opt; variable_space = Val(:log), verbose = 5, abstol = 1.0e-8, reltol = 1.0e-8))
sol = solve(prob, opt; variable_space = Val(:linear), verbose = 5, abstol = 1.0e-8, reltol = 1.0e-8)

function numpH(p, Vb)
    idx = findfirst(x -> first(x) === :Vb, p)
    if idx !== nothing
        p[idx] = :Vb => Vb
    end
    prob = EquilibriumProblem(A, μ, _n0(p)[indices], p = p)
    sol = solve(prob, opt; variable_space = Val(:linear), verbose = 0, abstol = 1.0e-8, reltol = 1.0e-8)
    pp = NamedTuple(p)
    Va, Vb = pp.Va, pp.Vb
    pHsol = sol.u[5] > sol.u[1] ?
        -log10(max(sol.u[5], pp.ϵ) / (Va + Vb)) :
        14 + log10(max(sol.u[1], pp.ϵ) / (Va + Vb))
    # pHsol = -log10(max(sol.u[5], pp.ϵ)/(Va+Vb))
    if !SciMLBase.successful_retcode(sol.retcode)
        println("Vb = ", Vb, ", pH = ", pHsol, " → success = ", SciMLBase.successful_retcode(sol.retcode), " → retcode = ", sol.retcode)
        println("  sol = ", sol.u)
        println("  sol.objective = ", sol.objective)
        # solve(prob, Val(:log), opt, verbose=5, abstol=1e-12, reltol=1e-12)
    end
    return pHsol
end

function anapH(p, Vb)
    pp = NamedTuple(p)
    ca, Va, cb, Vb = pp.ca, pp.Va, pp.cb, Vb
    V = Va + Vb
    Vbeq = ca * Va / cb
    return Vb < Vbeq ? pKa + log10(cb * Vb / (ca * Va - cb * Vb)) : 14 + log10((cb * Vb - ca * Va) / V)
end

pp = NamedTuple(p)
ca, Va, cb = pp.ca, pp.Va, pp.cb
Vbeq = ca * Va / cb

plot(xlabel = "Vᵇ", ylabel = "pH", legend = :none)
lVb = range(0.01Vbeq, 3Vbeq, 10000)
plot!(lVb, anapH.(Ref(p), lVb))
pH0 = (-log10(ca) + pKa) / 2
pHeq = 0.5 * (pKa + 14 + log10(ca * Va / (Va + Vbeq)))
scatter!([0, Vbeq], [pH0, pHeq])
lVb = range(0Vbeq, 3Vbeq, 200)
pHVb = numpH.(Ref(p), lVb)
scatter!(lVb, pHVb)


# function obj_ni(n, p)
#     pp = NamedTuple(p)
#     i = pp.i
#     return -n[i]
# end

# function cons!(res, n, p)
#     pp = NamedTuple(p)
#     res .= pp.A*n .- pp.b
# end

# u0 = _n0(p)[indices]
# b = A*u0

# optf = OptimizationFunction(
#         obj_ni,
#         Optimization.AutoForwardDiff();
#         cons = cons!
#     )

# prob = OptimizationProblem(optf, u0, [:A => A, :b => b, :i => 3];
#                            lcons = zeros(4), ucons = zeros(4),
#                            cons = cons!, lb = zeros(7), ub = fill(Inf, 7))
# sol = solve(prob, IpoptOptimizer(), print_level=0, abstol=1e-14, reltol=1e-18)
