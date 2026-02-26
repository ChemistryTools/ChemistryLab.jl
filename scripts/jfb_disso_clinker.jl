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
substances = build_species_from_database(df_substances)
init_species = split("C3S C2S Gp Anh Portlandite Jennite H2O@")
species = speciation(substances, init_species;
               aggregate_state=[AS_AQUEOUS],
               exclude_species=split("H2@ O2@ H2S@ HS- S2O3-2 SO3-2 HSO3-"),
               include_species=init_species)

cs = ChemicalSystem(species, CEMDATA_PRIMARIES)

state = ChemicalState(cs)
compo = ["C3S" => 67.8/100, "C2S" => 16.6/100, "Gp" => 2.8/100]
c =  sum(last.(compo))
wc = 0.4
w = wc*c
mtot = c+w
for x in compo set_quantity!(state, x.first, x.second/mtot*u"kg") end
set_quantity!(state, "H2O@", w/mtot*u"kg")
V = volume(state)
set_quantity!(state, "H+", 1.e-7*u"mol/L"*V.liquid)
set_quantity!(state, "OH-", 1.e-7*u"mol/L"*V.liquid)

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

solver = EquilibriumSolver(
    cs,
    DiluteSolutionModel(),
    opt;
    vartype = Val(:linear),
    abstol  = 1e-10,
    reltol  = 1e-10,
)
state_eq = solve(solver, state)
