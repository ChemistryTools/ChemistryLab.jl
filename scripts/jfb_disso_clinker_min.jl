using Pkg
Pkg.activate(@__DIR__)

using Revise
using ChemistryLab
using Optimization, OptimizationIpopt
using DynamicQuantities
using Optimization, OptimizationMOI, OptimizationOptimJL, OptimizationIpopt

substances = build_species("data/cemdata18-thermofun.json")
input_species = split("C3S C2S C3A C4AF Gp Anh Portlandite Jennite H2O@")
species = speciation(substances, input_species; aggregate_state = [AS_AQUEOUS])

cs = ChemicalSystem(species, CEMDATA_PRIMARIES)

# ── Initial state ──────────────────────────────────────────────────────────────
state = ChemicalState(cs)
compo = ["C3S" => 67.8 / 100, "C2S" => 16.6 / 100, "C3A" => 4 / 100, "C4AF" => 7.2 / 100, "Gp" => 2.8 / 100]
for x in compo
    set_quantity!(state, x.first, x.second * u"kg")
end
w = 0.4 * sum(last.(compo))
set_quantity!(state, "H2O@", w * u"kg")
V = volume(state)
set_quantity!(state, "H+", 1.0e-7u"mol/L" * V.liquid)
set_quantity!(state, "OH-", 1.0e-7u"mol/L" * V.liquid)
rescale!(state, 1.0u"kg")

# ── Equilibrate ─────────────────────────────────────────────────────────────────
state_eq = equilibrate(state; variable_space = Val(:linear), verbose = 5)

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
p = ChemistryLab._build_params(state; ϵ = 1.0e-16)
Gini = μ(ustrip.(state.n), p) ⋅ ustrip.(state.n)
Gfin = μ(ustrip.(state_eq.n), p) ⋅ ustrip.(state_eq.n)


values = [6.3267e+0, 5.2875e-6, 4.7599e-11, 6.0853e-10, 1.0e-16, 2.4934e-10, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16, 3.1801e-1, 1.0e-16, 6.5027e-13, 1.4563e-6, 4.3002e-1, 2.1376e-1, 1.0e-16, 1.0e-16, 6.9034e-8, 6.8628e-6, 7.9351e-3, 1.1895e-15, 6.9034e-8, 1.337e-3, 1.7474e-16, 1.0e-16, 8.0774e-16, 1.7135e-4, 1.0e-16, 1.0e-16, 1.0e-16, 2.5682e-7, 1.0e-16, 1.0e-16, 1.0e-16, 8.1095e-4, 5.8068e-9, 2.8552e+0, 1.0e-16, 1.0e-16, 1.0e-16, 1.0e-16, 3.534e+0, 1.1724e-1, 4.806e-16]
state_rkt = copy(state)
state_rkt.n .= values * u"mol"
ChemistryLab._update_derived!(state_rkt)
Grkt = μ(ustrip.(state_rkt.n), p) ⋅ ustrip.(state_rkt.n)
