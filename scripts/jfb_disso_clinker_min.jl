using Pkg
Pkg.activate(@__DIR__)

using Revise
using ChemistryLab
using DynamicQuantities
using Optimization, OptimizationMOI, OptimizationOptimJL, OptimizationIpopt

substances = build_species("data/cemdata18-thermofun.json")
input_species = split("C3S C2S C3A C4AF Gp Anh Portlandite Jennite H2O@")
species = speciation(substances, input_species; aggregate_state=[AS_AQUEOUS])

cs = ChemicalSystem(species, CEMDATA_PRIMARIES)

# ── Initial state ──────────────────────────────────────────────────────────────
state = ChemicalState(cs)
compo = ["C3S" => 67.8/100, "C2S" => 16.6/100, "C3A" => 4/100, "C4AF" => 7.2/100, "Gp" => 2.8/100]
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
state_eq = equilibrate(state; variable_space=Val(:linear), verbose=5)

display(state_eq)

# # ── Equilibrate with explicit options ───────────────────────────────────────────
# state_eq = equilibrate(state;
#     variable_space = Val(:linear),
#     abstol  = 1e-10,
#     reltol  = 1e-10,
# )

# display(state_eq)

# # ── Equilibrate with explicit solver ───────────────────────────────────────────
# opt = IpoptOptimizer(
#     acceptable_tol        = 1e-12,
#     dual_inf_tol          = 1e-12,
#     acceptable_iter       = 100,
#     constr_viol_tol       = 1e-12,
#     warm_start_init_point = "no",
# )

# state_eq = equilibrate(state;
#     solver  = opt,
#     variable_space = Val(:linear),
#     abstol  = 1e-10,
#     reltol  = 1e-10,
# )

# display(state_eq)

using LinearAlgebra

μ = build_potentials(cs, DiluteSolutionModel())
p  = ChemistryLab._build_params(state; ϵ=1.e-16)
Gini = μ(ustrip.(state.n), p) ⋅ ustrip.(state.n)
Gfin = μ(ustrip.(state_eq.n), p) ⋅ ustrip.(state_eq.n)


values = [6.3267e+00, 5.2875e-06, 4.7599e-11, 6.0853e-10, 1.0000e-16, 2.4934e-10, 1.0000e-16, 1.0000e-16, 1.0000e-16, 1.0000e-16, 1.0000e-16, 1.0000e-16, 1.0000e-16, 1.0000e-16, 3.1801e-01, 1.0000e-16, 6.5027e-13, 1.4563e-06, 4.3002e-01, 2.1376e-01, 1.0000e-16, 1.0000e-16, 6.9034e-08, 6.8628e-06, 7.9351e-03, 1.1895e-15, 6.9034e-08, 1.3370e-03, 1.7474e-16, 1.0000e-16, 8.0774e-16, 1.7135e-04, 1.0000e-16, 1.0000e-16, 1.0000e-16, 2.5682e-07, 1.0000e-16, 1.0000e-16, 1.0000e-16, 8.1095e-04, 5.8068e-09, 2.8552e+00, 1.0000e-16, 1.0000e-16, 1.0000e-16, 1.0000e-16, 3.5340e+00, 1.1724e-01, 4.8060e-16]
state_rkt = copy(state)
state_rkt.n .= values*u"mol"
ChemistryLab._update_derived!(state_rkt)
Grkt = μ(ustrip.(state_rkt.n), p) ⋅ ustrip.(state_rkt.n)
