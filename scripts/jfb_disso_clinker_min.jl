using Pkg
Pkg.activate(@__DIR__)

using Revise
using ChemistryLab
using DynamicQuantities
using Optimization, OptimizationMOI, OptimizationOptimJL, OptimizationIpopt

df_elements, df_substances, df_reactions = read_thermofun_database("data/cemdata18-thermofun.json")
substances = build_species_from_database(df_substances)
init_species = split("C3S C2S Gp Anh Portlandite Jennite H2O@")
species = speciation(substances, init_species;
               aggregate_state=[AS_AQUEOUS],
               exclude_species=split("H2@ O2@ H2S@ HS- S2O3-2 SO3-2 HSO3-"),
               include_species=init_species)

cs = ChemicalSystem(species, CEMDATA_PRIMARIES)

# ── Initial state ──────────────────────────────────────────────────────────────
state = ChemicalState(cs)
compo = ["C3S" => 67.8/100, "C2S" => 16.6/100, "Gp" => 2.8/100]
c     = sum(last.(compo))
wc    = 0.4
w     = wc * c
mtot  = c + w
for x in compo
    set_quantity!(state, x.first, x.second / mtot * u"kg")
end
set_quantity!(state, "H2O@", w / mtot * u"kg")
V = volume(state)
set_quantity!(state, "H+",  1e-7u"mol/L" * V.liquid)
set_quantity!(state, "OH-", 1e-7u"mol/L" * V.liquid)

# ── Equilibrate ─────────────────────────────────────────────────────────────────
state_eq = equilibrate(state)

display(state_eq)

# ── Equilibrate with explicit options ───────────────────────────────────────────
state_eq = equilibrate(state;
    vartype = Val(:linear),
    abstol  = 1e-10,
    reltol  = 1e-10,
)

display(state_eq)

# ── Equilibrate with explicit solver ───────────────────────────────────────────
opt = IpoptOptimizer(
    acceptable_tol        = 1e-12,
    dual_inf_tol          = 1e-12,
    acceptable_iter       = 100,
    constr_viol_tol       = 1e-12,
    warm_start_init_point = "no",
)

state_eq = equilibrate(state;
    solver  = opt,
    vartype = Val(:linear),
    abstol  = 1e-10,
    reltol  = 1e-10,
)

display(state_eq)
