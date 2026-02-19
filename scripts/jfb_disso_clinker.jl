using Pkg
Pkg.activate(@__DIR__)
# Pkg.add(path="../ChemistryLab.jl")

using Revise
using ChemistryLab
using DynamicQuantities
using Optimization, OptimizationMOI, OptimizationOptimJL, OptimizationIpopt
using Plots
using SparseArrays
using LinearAlgebra

df_elements, df_substances, df_reactions = read_thermofun_database("data/cemdata18-thermofun.json")
dict_species = build_species_from_database(df_substances)

df_selection = get_compatible_species(df_substances, split("C3S C2S Gp Anh Portlandite Jennite H2O@");
               aggregate_states=[AS_AQUEOUS], exclude_species=split("H2@ O2@ H2S@ HS- S2O3-2 SO3-2 HSO3-"), union=true)
d_species = build_species_from_database(df_selection)
primaries = collect(skipmissing(get.(Ref(d_species), CEMDATA_PRIMARIES, missing)))

species = collect(values(d_species))

sort!(species, by=x->(x.aggregate_state, x.class, x.symbol))

SM = StoichMatrix(species, primaries) ; pprint(SM)

μ = potentials_dilute_ideal(species)

ρ(s) = ustrip(s.M)/s.V⁰()

calc_ρᶜ(compo) = sum(x -> x.second, compo) / sum(x -> x.second/ρ(d_species[x.first]), compo)

function initial_content(cem_compo, cem_formul, ρʷ, ρᵍ)
    ρᶜ = calc_ρᶜ(cem_compo)
    dict_cem_formul = Dict(cem_formul)
    wc = get(dict_cem_formul, :wc, 0.)
    gc = get(dict_cem_formul, :gc, 0.)
    ϕᵃⁱʳ = get(dict_cem_formul, :ϕᵃⁱʳ, 0.)
    N = Pair{String,Number}[elem.first => elem.second * ρᶜ * (1-ϕᵃⁱʳ)/ustrip(d_species[elem.first].M)/(1+wc*ρᶜ/ρʷ+gc*ρᶜ/ρᵍ) for elem ∈ cem_compo]
    NGp = get(Dict(N), "Gp", 0)
    push!(N, "H2O@" => wc * ρᶜ * (1-ϕᵃⁱʳ)/ustrip(d_species["H2O@"].M)/(1+wc*ρᶜ/ρʷ+gc*ρᶜ/ρᵍ) - 3/2 * NGp)
    # push!(N, :ϕ => ϕᵃⁱʳ + wc*ρᶜ/ρʷ * (1-ϕᵃⁱʳ)/(1+wc*ρᶜ/ρʷ+gc*ρᶜ/ρᵍ))
    ϕ = ϕᵃⁱʳ + wc*ρᶜ/ρʷ * (1-ϕᵃⁱʳ)/(1+wc*ρᶜ/ρʷ+gc*ρᶜ/ρᵍ)
    return N, ϕ
end

# compo = ["C3S" => 67.8/100, "C2S" => 16.6/100, "C3A" => 4/100, "C4AF" => 7.2/100, "Gp" => 2.8/100]
compo = ["C3S" => 67.8/100, "C2S" => 16.6/100, "Gp" => 2.8/100]
ρʷ = 1000.0 ; ρᵍ = 2700.0 ;
N, ϕ = initial_content(compo, [:wc => 0.4, :gc => 0., :ϕᵃⁱʳ => 0.], ρʷ, ρᵍ)
N = Dict(N)
Vinit = ϕ*1000
n₀ = [get(N, k, k=="OH-" || k=="H+" ? 1e-7*Vinit : 0.) for k in symbol.(species)]

T = 298.15
p = [:T => T, :ϵ => 1.e-16, :ΔₐG⁰ => [s.ΔₐG⁰(T = T) for s in species]]

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
prob = EquilibriumProblem(Float64.(SM.A), μ, n₀, p=p)
sol = solve(prob, opt, Val(:linear); verbose=5, abstol=1e-10, reltol=1e-10)

initial = Dict(s => n₀[i] for (i,s) in enumerate(symbol.(species)))
final = Dict(s => sol[i] for (i,s) in enumerate(symbol.(species)))
V = final["H2O@"]*ustrip(d_species["H2O@"].M)
h = final["H+"]/V
ω = final["OH-"]/V

Gini = n₀ ⋅ μ(n₀, p)
Gfin = sol ⋅ μ(sol, p)

reaktoro = Dict{String, Float64}(
    "H+"          => 1.0117e-10,
    "SiO2@"       => 4.8240e-7,
    "Portlandite" => 5.9390e+3,
    "OH-"         => 5.4832e+0,
    "SO4-2"       => 5.3145e-1,
    "CaSiO3@"     => 8.7804e-3,
    "Anh"         => 1.9468e-16,
    "SiO3-2"      => 1.8941e-5,
    "Jennite"     => 5.4601e+3,
    "Ca+2"        => 2.7433e+0,
    "Gp"          => 2.2400e+2,
    "HSO4-"       => 2.2185e-11,
    "C3S"         => 1.0000e-16,
    "C2S"         => 1.0000e-16,
    "Ca(SO4)@"    => 1.2349e+0,
    "HSiO3-"      => 1.7391e-4,
    "CaOH+"       => 1.0596e+0,
    "H2O@"        => 1.3078e+4,
    "Si4O10-4"    => 1.0000e-16,
    "Ca(HSiO3)+"  => 3.2098e-5
)
sol_rkt = [reaktoro[s] for s in symbol.(species)]
Grkt = sol_rkt ⋅ μ(sol_rkt, p)

prob2 = EquilibriumProblem(Float64.(SM.A), μ, sol_rkt, p=p)
sol2 = solve(prob, opt, Val(:linear); verbose=5, abstol=1e-10, reltol=1e-10)
final2 = Dict(s => sol2[i] for (i,s) in enumerate(symbol.(species)))
