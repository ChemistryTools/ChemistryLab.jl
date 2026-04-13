# =============================================================================
# lavergne2018_cement_calorimetry.jl
#
# Simulation du dégagement de chaleur d'une pâte de CEM I en calorimètre
# semi-adiabatique par la méthode de Lavergne et al. (2018).
#
# Cinétique d'hydratation : modèle de Parrot & Killoh (1984)
#   via ChemistryLab.parrot_killoh (KineticFunc),
#   température : correction d'Arrhenius par phase,
#   eau disponible : limite maximale d'hydratation α_max (Powers 1948).
#
# Calorimètre : semi-adiabatique, perte de chaleur quadratique
#   φ(ΔT) = a·ΔT + b·ΔT²
#   (ChemistryLab.SemiAdiabaticCalorimeter)
#
# Références :
#   Lavergne F., Ben Fraj A., Bayane I., Barthélémy J.-F. (2018).
#   Estimating the mechanical properties of hydrating blended cementitious
#   materials: an investigation based on micromechanics.
#   Cem. Concr. Res. 104, 37-60.
#   https://doi.org/10.1016/j.cemconres.2017.11.007
#
#   Parrot L.J. & Killoh D.C. (1984).
#   Prediction of cement hydration.
#   Br. Ceram. Proc. 35, 41-53.
#
# Usage :
#   julia --project scripts/lavergne2018_cement_calorimetry.jl
#   ou depuis le REPL : include("scripts/lavergne2018_cement_calorimetry.jl")
# =============================================================================

using Pkg
Pkg.activate(@__DIR__)

using ChemistryLab
using OrdinaryDiffEq
using DynamicQuantities
using OrderedCollections
using Printf

# ── 1. ChemicalSystem depuis la base CEMDATA18 ────────────────────────────────

const DATA_FILE = joinpath(pkgdir(ChemistryLab), "data", "cemdata18-thermofun.json")

substances = build_species(DATA_FILE)

input_species = split(
    "C3S C2S C3A C4AF " *
        "Portlandite Jennite ettringite monosulphate12 C3AH6 C3FH6 " *
        "H2O@",
)

species = speciation(substances, input_species; aggregate_state = [AS_AQUEOUS])
cs = ChemicalSystem(species, CEMDATA_PRIMARIES)

@info "Système chimique : $(length(cs.species)) espèces, $(length(cs.reactions)) réactions"

# ── 2. Composition de la pâte de ciment ──────────────────────────────────────
#
# CEM I 52,5 type, d'après Lavergne et al. (2018)
# Fractions massiques dans le clinker

const PHASE_MASS_FRAC = (C3S = 0.619, C2S = 0.165, C3A = 0.08, C4AF = 0.087)

# Rapport eau/ciment
const WC = 0.4

# Degré maximal d'hydratation (loi de Powers 1948 : α_max ≤ w/c / 0.42)
const ALPHA_MAX = min(1.0, WC / 0.42)

state0 = ChemicalState(cs)
for (name, frac) in pairs(PHASE_MASS_FRAC)
    set_quantity!(state0, string(name), frac * u"kg")
end
set_quantity!(state0, "H2O@", WC * u"kg")

# ── 3. Modèles cinétiques de Parrot & Killoh ──────────────────────────────────
#
# Paramètres K₁, N₁, K₂, N₂, K₃, N₃, B, Ea — Parrot & Killoh (1984)
# Ea — Schindler & Folliard (2005) / van Breugel (1991) [J/mol]

const PK_L2018_C3S = (
    K₁ = 1.5u"1/d", N₁ = 3.3, K₂ = 0.018u"1/d", N₂ = 2.5,
    K₃ = 0.0024u"1/d", N₃ = 4.0, B = 0.5,
    Ea = 41_570.0u"J/mol", T_ref = 293.15u"K",
)
const PK_L2018_C2S = (
    K₁ = 0.95u"1/d", N₁ = 0.5, K₂ = 0.0005u"1/d", N₂ = 2.5,
    K₃ = 0.0024u"1/d", N₃ = 4.0, B = 0.2,
    Ea = 43_670.0u"J/mol", T_ref = 293.15u"K",
)
const PK_L2018_C3A = (
    K₁ = 0.082u"1/d", N₁ = 0.87, K₂ = 0.00024u"1/d", N₂ = 2.0,
    K₃ = 0.0024u"1/d", N₃ = 4.0, B = 0.04,
    Ea = 54_040.0u"J/mol", T_ref = 293.15u"K",
)
const PK_L2018_C4AF = (
    K₁ = 0.165u"1/d", N₁ = 3.7, K₂ = 0.0015u"1/d", N₂ = 2.5,
    K₃ = 0.0024u"1/d", N₃ = 4.0, B = 0.5,
    Ea = 34_420.0u"J/mol", T_ref = 293.15u"K",
)

pk_C3S = parrot_killoh(PK_L2018_C3S, "C3S"; α_max = ALPHA_MAX)
pk_C2S = parrot_killoh(PK_L2018_C2S, "C2S"; α_max = ALPHA_MAX)
pk_C3A = parrot_killoh(PK_L2018_C3A, "C3A"; α_max = ALPHA_MAX)
pk_C4AF = parrot_killoh(PK_L2018_C4AF, "C4AF"; α_max = ALPHA_MAX)

# ── 4. Liste des réactions cinétiques ─────────────────────────────────────────
#
# Réactions d'hydratation simplifiées (CEMDATA18 / Lothenbach & Winnefeld 2006),
# coefficients stœchiométriques calculés pour Jennite = Ca₉Si₆O₁₈(OH)₆·8H₂O.
#
#   C₃S  + 3.33 H₂O  →  0.167 Jennite + 1.5 Portlandite   (Ca:Si = 1.5 dans Jennite)
#   C₂S  + 2.33 H₂O  →  0.167 Jennite + 0.5 Portlandite
#   C₃A  + 6    H₂O  →  C₃AH₆
#   C₄AF + 6    H₂O  →  0.5 C₃AH₆ + 0.5 C₃FH₆ + Ca(OH)₂  (approx.)
#
# La chaleur est calculée automatiquement depuis les ΔₐH⁰ des espèces CEMDATA18
# (via complete_thermo_functions! / _reaction_enthalpy).

sp(name) = cs[name]

rxn_C3S = Reaction(
    OrderedDict(sp("C3S") => 1.0, sp("H2O@") => 3.33),
    OrderedDict(sp("Jennite") => 0.167, sp("Portlandite") => 1.5);
    symbol = "C₃S hydration",
)
rxn_C3S[:rate] = pk_C3S

rxn_C2S = Reaction(
    OrderedDict(sp("C2S") => 1.0, sp("H2O@") => 2.33),
    OrderedDict(sp("Jennite") => 0.167, sp("Portlandite") => 0.5);
    symbol = "C₂S hydration",
)
rxn_C2S[:rate] = pk_C2S

rxn_C3A = Reaction(
    OrderedDict(sp("C3A") => 1.0, sp("H2O@") => 6.0),
    OrderedDict(sp("C3AH6") => 1.0);
    symbol = "C₃A hydration",
)
rxn_C3A[:rate] = pk_C3A

rxn_C4AF = Reaction(
    OrderedDict(sp("C4AF") => 1.0, sp("H2O@") => 6.0),
    OrderedDict(sp("C3AH6") => 0.5, sp("C3FH6") => 0.5, sp("Portlandite") => 1.0);
    symbol = "C₄AF hydration",
)
rxn_C4AF[:rate] = pk_C4AF

kinetic_reactions = [rxn_C3S, rxn_C2S, rxn_C3A, rxn_C4AF]

# ── 6. Calorimètre semi-adiabatique ──────────────────────────────────────────
#
# Dispositif de type Langavant (norme NF EN 196-9).
# Perte de chaleur quadratique φ(ΔT) = a·ΔT + b·ΔT² (Lavergne et al. 2018).
# Capacité thermique Cp [J/K] pour 1 kg ciment + WC kg eau + vase :
#   cp_ciment ≈ 800 J/(kg·K), cp_eau = 4186 J/(kg·K), cp_vase ≈ 900 J/(kg·K)

const T0_K = 20.0 + 273.15    # température initiale [K]
const T_ENV_K = 20.0 + 273.15  # température ambiante [K]
const CP = 1.0 * 800.0 + WC * 4186.0 + 1.0 * 900.0   # ≈ 3449 J/K

const A_CAL = 0.3   # coefficient linéaire  [W/K]
const B_CAL = 0.003  # coefficient quadratique [W/K²]

cal = SemiAdiabaticCalorimeter(;
    Cp = CP * u"J/K",
    T_env = T_ENV_K * u"K",
    heat_loss = ΔT -> A_CAL * ΔT + B_CAL * ΔT^2,
    T0 = T0_K * u"K",
)

# ── 7. Problème cinétique ─────────────────────────────────────────────────────

const TSPAN = (0.0, 7.0 * 86400.0)

kp = KineticsProblem(
    cs,
    kinetic_reactions,
    state0,
    TSPAN;
    equilibrium_solver = nothing,
)

# ── 8. Intégration ────────────────────────────────────────────────────────────

@info "Résolution en cours (7 jours, Rodas5P)..."
ks = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-9)
sol = integrate(kp, ks; calorimeter = cal)
@info "Résolution terminée : $(length(sol.t)) pas acceptés."

# ── 9. Post-traitement ────────────────────────────────────────────────────────

t_h = sol.t ./ 3600.0   # secondes → heures

t_T, T_K_vec = temperature_profile(sol, cal)
t_Q, Q_J_vec = cumulative_heat(sol, cal)
T_°C_vec = T_K_vec .- 273.15
Q_kJ_vec = Q_J_vec ./ 1000.0

_, qdot_W = heat_flow(sol, cal)   # [W]

# Degrés d'hydratation par phase
n0_kin = [sol.prob.p.n_initial_full[i] for i in kp.idx_kin_unique]
n_kin = [[u[i] for u in sol.u] for i in eachindex(n0_kin)]

# idx_kin_unique liste les espèces cinétiques dans l'ordre d'apparition
# La position de chaque phase clinker dépend de cs.species; on la récupère
# en comparant avec kp.idx_kin_unique.
function phase_alpha(cs, kp, sol, n0_kin, n_kin, name)
    sp_idx = findfirst(sp -> ChemistryLab.phreeqc(ChemistryLab.formula(sp)) == name, cs.species)
    pos = findfirst(==(sp_idx), kp.idx_kin_unique)
    isnothing(pos) && return fill(NaN, length(sol.t))
    return 1.0 .- n_kin[pos] ./ n0_kin[pos]
end

α_C3S_vec = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C3S")
α_C2S_vec = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C2S")
α_C3A_vec = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C3A")
α_C4AF_vec = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C4AF")

# Degré d'hydratation moyen pondéré par la fraction massique
α_mean = (
    PHASE_MASS_FRAC.C3S .* α_C3S_vec .+
        PHASE_MASS_FRAC.C2S .* α_C2S_vec .+
        PHASE_MASS_FRAC.C3A .* α_C3A_vec .+
        PHASE_MASS_FRAC.C4AF .* α_C4AF_vec
)

# ── Résumé ────────────────────────────────────────────────────────────────────
println()
println("╔══════════════════════════════════════════════╗")
println("║    Résultats calorimétrie semi-adiabatique   ║")
println("║    CEM I — Lavergne et al. (2018)            ║")
println("╠══════════════════════════════════════════════╣")
println("║  Temps       = 7 jours")
@printf "║  T final     = %.2f °C\n" T_°C_vec[end]
@printf "║  ΔT max      = %.2f °C\n" maximum(T_°C_vec) - (T0_K - 273.15)
@printf "║  α(C3S)      = %.4f\n"   α_C3S_vec[end]
@printf "║  α(C2S)      = %.4f\n"   α_C2S_vec[end]
@printf "║  α(C3A)      = %.4f\n"   α_C3A_vec[end]
@printf "║  α(C4AF)     = %.4f\n"   α_C4AF_vec[end]
@printf "║  α moyen     = %.4f\n"   α_mean[end]
@printf "║  Q total     = %.1f kJ/kg_ciment\n" Q_kJ_vec[end]
println("╚══════════════════════════════════════════════╝")

# ── 10. Tracés ─────────────────────────────────────────────────────────────────

using Plots
gr()

p1 = plot(
    t_T ./ 3600, T_°C_vec;
    xlabel = "Temps [h]", ylabel = "Température [°C]",
    title = "Profil de température — calorimètre semi-adiabatique",
    label = "T(t)", lw = 2, color = :red,
)
hline!(p1, [T0_K - 273.15]; linestyle = :dash, color = :gray, label = "T₀ = T_env")

p2 = plot(
    t_h, [α_C3S_vec α_C2S_vec α_C3A_vec α_C4AF_vec α_mean];
    xlabel = "Temps [h]", ylabel = "Degré d'hydratation α",
    title = "Hydratation des phases clinker (modèle Parrot-Killoh)",
    label = ["C₃S" "C₂S" "C₃A" "C₄AF" "ᾱ moyen"],
    lw = 2, ls = [:solid :dash :dot :dashdot :solid],
)
hline!(p2, [ALPHA_MAX]; linestyle = :dash, color = :black, label = "α_max")

p3 = plot(
    t_T ./ 3600, qdot_W ./ 1000;
    xlabel = "Temps [h]", ylabel = "Flux de chaleur [kW/kg_ciment]",
    title = "Flux thermique instantané",
    label = "q̇(t)", lw = 2, color = :orange,
)

p4 = plot(
    t_Q ./ 3600, Q_kJ_vec;
    xlabel = "Temps [h]", ylabel = "Q [kJ/kg_ciment]",
    title = "Chaleur cumulée",
    label = "Q(t)", lw = 2, color = :purple,
)

display(
    plot(
        p1, p2, p3, p4; layout = (2, 2), size = (1200, 800),
        plot_title = "Lavergne et al. (2018) — CEM I, E/C = $WC",
    )
)
