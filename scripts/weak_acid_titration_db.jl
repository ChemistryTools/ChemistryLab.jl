using Pkg
Pkg.activate(@__DIR__)
# Pkg.add(path="../ChemistryLab.jl")

using Revise
using ChemistryLab
using DynamicQuantities
using Optimization, OptimizationMOI, OptimizationOptimJL, OptimizationIpopt
using Plots
using SparseArrays

df_elements_orga, df_substances_orga, df_reactions_orga = read_thermofun_database("data/slop98-organic-thermofun.json")
dict_species_orga = Dict(symbol(s) => s for s in build_species_from_database(df_substances_orga))
df_elements_inorga, df_substances_inorga, df_reactions_inorga = read_thermofun_database("data/slop98-inorganic-thermofun.json")
dict_species_inorga = Dict(symbol(s) => s for s in build_species_from_database(df_substances_inorga))
dict_species = merge(dict_species_orga, dict_species_inorga)

species = [dict_species[s] for s in split("H2O@ OH- H+ AceH@ Ace- NaOH@ Na+")]
H₂O, OH⁻, H⁺, AH, A⁻, BOH, B⁺ = species
primaries = [H₂O, H⁺, A⁻, B⁺]
SM = StoichMatrix(species, primaries) ; pprint(SM)

μ = potentials_dilute_ideal(species)

function n0(p)
    p = NamedTuple(p)
    ϵ = p.ϵ
    order = p.order
    ca, Va, cb, Vb = p.ca, p.Va, p.cb, p.Vb
    V = Va + Vb

    n = Vector{Float64}(undef, length(species))

    n[order["H2O@"]] = 55.5*V
    n[order["H+"]] = 1.e-7V
    n[order["OH-"]] = 1.e-7V

    n[order["AceH@"]] = max(ca*Va, ϵ)
    n[order["Ace-"]] = ϵ
    n[order["NaOH@"]] = max(cb*Vb, ϵ)
    n[order["Na+"]] = ϵ

    return n
end

function ub(p)
    p = NamedTuple(p)
    ϵ = p.ϵ
    order = p.order
    ca, Va, cb, Vb = p.ca, p.Va, p.cb, p.Vb
    V = Va + Vb

    v = Vector{Float64}(undef, length(species))

    v[order["H2O@"]] = 2*55.5V
    v[order["H+"]] = 1.e-1V
    v[order["OH-"]] = 1.e-1V

    v[order["AceH@"]] = max(ca*Va, ϵ)
    v[order["Ace-"]] = max(ca*Va, ϵ)
    v[order["NaOH@"]] = max(cb*Vb, ϵ)
    v[order["Na+"]] = max(cb*Vb, ϵ)

    return v
end

R = ustrip(Constants.R)
T = 298.15
RT = R*T

# Data acide acétique CH₃COOH et soude NaOH
Ka = exp(-(A⁻.ΔₐG⁰(T = T)+H⁺.ΔₐG⁰(T = T)-AH.ΔₐG⁰(T = T))/RT)
pKa = -log10(Ka)
Kb = exp(-(B⁺.ΔₐG⁰(T = T)+OH⁻.ΔₐG⁰(T = T)-BOH.ΔₐG⁰(T = T))/RT)
pKb = -log10(Kb)

order_species =  Dict(symbol(species[i]) => i for i in eachindex(species))
p = [:order => order_species,
     :T => T,
     :ΔₐG⁰ => [s.ΔₐG⁰(T = T) for s in species],
     :idxsolvent => findfirst(x->x.class == SC_AQSOLVENT, species),
     :ϵ => 1.e-16, :ca => 0.1, :Va => 0.1, :cb => 0.1, :Vb => 0.]

prob = EquilibriumProblem(SM.A, μ, n0(p), ub=ub(p), p=p)
opt = IpoptOptimizer(
    acceptable_tol = 1e-12,
    dual_inf_tol = 1e-12,
    acceptable_iter = 100,
    constr_viol_tol = 1e-12,
    # compl_inf_tol = 1e-4,
    # mu_strategy = "adaptive",
    warm_start_init_point = "no",
    # expect_infeasible_problem = "yes",
)
sol = exp.(solve(prob, opt, Val(:log), verbose=5, abstol=1e-10, reltol=1e-10))
sol = solve(prob, opt, Val(:linear), verbose=5, abstol=1e-10, reltol=1e-10)

idxOH⁻ = order_species["OH-"]
idxH⁺ = order_species["H+"]

function numpH(p, Vb)
    idx = findfirst(x -> first(x) === :Vb, p)
    if idx !== nothing p[idx] = :Vb => Vb end
    prob = EquilibriumProblem(SM.A, μ, n0(p), ub=ub(p), p=p)
    sol = solve(prob, opt, Val(:linear), verbose=0, abstol=1e-12, reltol=1e-12)
    pp = NamedTuple(p)
    Va, Vb = pp.Va, pp.Vb
    V = Va + Vb
    pHsol = sol.u[idxH⁺]>sol.u[idxOH⁻] ?
            -log10(max(sol.u[idxH⁺], pp.ϵ)/V) :
            14+log10(max(sol.u[idxOH⁻], pp.ϵ)/V)
    # pHsol = -log10(max(sol.u[idxH⁺], pp.ϵ)/V)
    if !SciMLBase.successful_retcode(sol.retcode)
        println("Vb = ", Vb,", pH = ", pHsol," → success = ", SciMLBase.successful_retcode(sol.retcode), " → retcode = ", sol.retcode)
        println("  sol = ", sol.u)
        println("  sol.objective = ", sol.objective)
    end
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
