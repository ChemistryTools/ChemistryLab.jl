using Pkg
Pkg.activate(@__DIR__)
# Pkg.add(path="../ChemistryLab.jl")

using Revise, ChemistryLab
using OptimizationMOI, OptimizationOptimJL, Ipopt
using Plots

R, T = 8.314, 273.15+25
RT = R*T

function μ(n, p)
    pp = NamedTuple(p)
    _n = max.(n, pp.ϵ)
    ΔGf₀ = pp.ΔGf₀
    ntot = sum(_n)
    nden = [ntot ; fill(n[1], length(n)-1)]
    return last.(ΔGf₀)/RT + log.(_n ./ nden) + [first(x) === :H2O ? 0. : log(55.5) for x in ΔGf₀]
end

function _n0(p)
    pp = NamedTuple(p)
    ϵ = pp.ϵ
    ca, Va, cb, Vb = pp.ca, pp.Va, pp.cb, pp.Vb
    V = Va+Vb

    n_H2O = 55.5*V
    n_H = 1.e-7V
    n_OH = 1.e-7V

    n_AH = max(ca*Va, ϵ)
    n_A = ϵ
    n_BOH = max(cb*Vb, ϵ)
    n_B = ϵ

    return [n_H2O, n_H, n_OH, n_AH, n_A, n_BOH, n_B]
end

A = [
    2 1 1 1 0 1 0
    1 0 1 0 0 1 0
    0 0 0 1 1 0 0
    0 0 0 0 0 1 1
]

function _ub(p)
    pp = NamedTuple(p)
    ϵ = pp.ϵ
    ca, Va, cb, Vb = pp.ca, pp.Va, pp.cb, pp.Vb
    V = Va+Vb

    return [2*55.5V, 1.e-1V, 1.e-1V, max(ca*Va, ϵ), max(ca*Va, ϵ), max(cb*Vb, ϵ), max(cb*Vb, ϵ)]
end

# Data acide acétique CH₃COOH et soude NaOH
ΔGf₀ = [
    :H2O => -237.18e3,
    :H => 0.e3,
    :OH => -157.27e3,
    :AH => -396.46e3,
    :A => -369.31e3,
    :BOH => -379.5e3,
    :B => -261.9e3
]
_ΔGf₀ = NamedTuple(ΔGf₀)
Ka = exp(-(_ΔGf₀[:A]+_ΔGf₀[:H]-_ΔGf₀[:AH])/RT)
pKa = -log10(Ka)
Kb = exp(-(_ΔGf₀[:B]+_ΔGf₀[:OH]-_ΔGf₀[:BOH])/RT)
pKb = -log10(Kb)

dep_comp = first.(ΔGf₀)
indep_comp = [:H, :O, :A, :B]

p = [:ΔGf₀ => ΔGf₀, :ϵ => 1.e-16, :ca => 0.1, :Va => 0.1, :cb => 0.1, :Vb => 0.]

prob = EquilibriumProblem(A, μ, _n0(p), ub=_ub(p), p=p, indep_components=indep_comp, dep_components=dep_comp)
sol = solve(prob, Ipopt.Optimizer(), print_level=0, tol=1e-16)

function numpH(p, Vb)
    idx = findfirst(x -> first(x) === :Vb, p)
    if idx !== nothing p[idx] = :Vb => Vb end
    prob = EquilibriumProblem(A, μ, _n0(p), ub=_ub(p), p=p, indep_components=indep_comp, dep_components=dep_comp)
    sol = solve(prob, Ipopt.Optimizer(), print_level=0, tol=1e-16)
    pp = NamedTuple(p)
    Va, Vb = pp.Va, pp.Vb
    pHsol = -log10(max(sol.u[2], pp.ϵ)/(Va+Vb))
    return pHsol
end

function anapH(p, Vb)
    pp = NamedTuple(p)
    ca, Va, cb, Vb = pp.ca, pp.Va, pp.cb, Vb
    V = Va+Vb
    Vbeq = ca*Va/cb
    return Vb<Vbeq ? pKa+log10(cb*Vb/(ca*Va-cb*Vb)) : 14+log10((cb*Vb-ca*Va)/V)
end

pp = NamedTuple(p)
ca, Va, cb = pp.ca, pp.Va, pp.cb
Vbeq = ca*Va/cb

plot(xlabel="Vᵇ", ylabel="pH", legend=:none)
lVb = range(0.01Vbeq, 3Vbeq, 10000)
plot!(lVb, anapH.(Ref(p), lVb))
pH0 = (-log10(ca)+pKa)/2
pHeq = 0.5*(pKa+14+log10(ca*Va/(Va+Vbeq)))
scatter!([0, Vbeq], [pH0, pHeq])
lVb = range(0Vbeq, 3Vbeq, 200)
pHVb = numpH.(Ref(p), lVb)
scatter!(lVb, pHVb)
