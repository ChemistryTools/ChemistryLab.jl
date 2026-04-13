# =============================================================================
# cement_clinker_kinetics.jl
#
# Simulation cinétique de l'hydratation du clinker OPC via ChemistryLab :
#   - ChemicalSystem construit depuis la base de données CEMDATA18
#   - KineticsProblem avec parrot_killoh (KineticFunc) pour les 4 phases clinker
#   - Suivi de la chaleur via SemiAdiabaticCalorimeter
#
# Ce script montre l'utilisation de la chaîne complète :
#   1. ChemicalSystem  → 2. ChemicalState  → 3. liste de Reaction annotées
#   → 4. KineticsProblem → 5. integrate → 6. post-traitement
#
# Cinétique : Parrot & Killoh (1984), correction Arrhenius de Schindler &
#   Folliard (2005).
# Calorimètre semi-adiabatique avec pertes quadratiques (Lavergne et al. 2018).
#
# Usage :
#   julia --project scripts/cement_clinker_kinetics.jl
#   ou depuis le REPL : include("scripts/cement_clinker_kinetics.jl")
# =============================================================================

using Pkg
Pkg.activate(@__DIR__)

using ChemistryLab
using OrdinaryDiffEq
using DynamicQuantities
using OrderedCollections
using Printf

# ── 1. ChemicalSystem depuis la base CEMDATA18 ────────────────────────────────
#
# On sélectionne :
#   • les 4 phases clinker (espèces cinétiques)
#   • les principaux produits d'hydratation
#   • l'eau liquide (solvant)
# L'extension CEMDATA_PRIMARIES donne les composantes indépendantes du système.

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

# ── 2. État initial ───────────────────────────────────────────────────────────
#
# Composition d'un CEM I 52,5 R typique (Lavergne et al. 2018) :
# fractions massiques dans le clinker. On travaille avec 1 kg de ciment.

const WC = 0.4         # rapport eau/ciment [-]
const COMPOSITION = (            # fractions massiques des phases clinker [kg/kg ciment]
    C3S = 0.619,
    C2S = 0.165,
    C3A = 0.08,
    C4AF = 0.087,
)

state0 = ChemicalState(cs)

for (name, frac) in pairs(COMPOSITION)
    set_quantity!(state0, string(name), frac * u"kg")
end
set_quantity!(state0, "H2O@", WC * u"kg")

# ── 3. Modèles cinétiques de Parrot & Killoh ──────────────────────────────────
#
# Degré maximal d'hydratation selon Powers (1948) : α_max ≤ w/c / 0.42
# On crée des modèles avec le α_max calculé plutôt que la valeur par défaut 1.0.

const α_max = min(1.0, WC / 0.42)

pk_C3S = parrot_killoh(PK_PARAMS_C3S, "C3S"; α_max)
pk_C2S = parrot_killoh(PK_PARAMS_C2S, "C2S"; α_max)
pk_C3A = parrot_killoh(PK_PARAMS_C3A, "C3A"; α_max)
pk_C4AF = parrot_killoh(PK_PARAMS_C4AF, "C4AF"; α_max)

# ── 4. Liste des réactions cinétiques ─────────────────────────────────────────
#
# Réactions d'hydratation (CEMDATA18 / Lothenbach & Winnefeld 2006) :
#   Jennite = Ca₉Si₆O₁₈(OH)₆·8H₂O  →  Ca:Si = 1.5
#
#   C₃S  + 3.33 H₂O  →  0.167 Jennite + 1.5  Portlandite
#   C₂S  + 2.33 H₂O  →  0.167 Jennite + 0.5  Portlandite
#   C₃A  + 6    H₂O  →  C₃AH₆
#   C₄AF + 6    H₂O  →  0.5 C₃AH₆ + 0.5 C₃FH₆ + Portlandite  (approx.)
#
# Chaque Reaction porte son taux cinétique (:rate) et son enthalpie (:heat_per_mol).
# La stœchiométrie des produits permet le suivi des phases hydratées.
#
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

# ── 5. Problème cinétique + calorimètre ──────────────────────────────────────
#
# tspan : 7 jours en secondes
# equilibrium_solver = nothing → pas de re-spéciation aqueuse (suffisant pour PK)

const TSPAN = (0.0, 7.0 * 86400.0)

kp = KineticsProblem(
    cs,
    kinetic_reactions,
    state0,
    TSPAN;
    equilibrium_solver = nothing,
)

# Paramètres calorimètre semi-adiabatique (Lavergne et al. 2018)
#   Cp [J/K] : 1 kg ciment + WC kg eau + vase Dewar
#   pertes quadratiques : φ(ΔT) = a·ΔT + b·ΔT²
cal = SemiAdiabaticCalorimeter(;
    Cp = (1.0 * 800.0 + WC * 4186.0 + 1.0 * 900.0) * u"J/K",   # ≈ 3449 J/K
    T_env = 293.15u"K",
    heat_loss = ΔT -> 0.3 * ΔT + 0.003 * ΔT^2,
    T0 = 293.15u"K",
)

# ── 6. Intégration ────────────────────────────────────────────────────────────

@info "Intégration en cours (7 jours, Rodas5P) …"
ks = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-9)
sol = integrate(kp, ks; calorimeter = cal)
@info "Terminé : $(length(sol.t)) pas acceptés."

# ── 7. Post-traitement ────────────────────────────────────────────────────────

t_h = sol.t ./ 3600.0   # [s] → [h]

# Température et chaleur depuis le calorimètre
t_T, T_K_vec = temperature_profile(sol, cal)
t_Q, Q_J_vec = cumulative_heat(sol, cal)
T_°C_vec = T_K_vec .- 273.15
Q_kJ_vec = Q_J_vec ./ 1000.0

# Degrés d'hydratation par phase
n0_kin = [sol.prob.p.n_initial_full[i] for i in kp.idx_kin_unique]
n_kin = [[u[i] for u in sol.u] for i in eachindex(n0_kin)]

function phase_alpha(cs, kp, sol, n0_kin, n_kin, name)
    sp_idx = findfirst(sp -> ChemistryLab.phreeqc(ChemistryLab.formula(sp)) == name, cs.species)
    pos = findfirst(==(sp_idx), kp.idx_kin_unique)
    isnothing(pos) && return fill(NaN, length(sol.t))
    return 1.0 .- n_kin[pos] ./ n0_kin[pos]
end

α_C3S = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C3S")
α_C2S = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C2S")
α_C3A = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C3A")
α_C4AF = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C4AF")

# Degré d'hydratation moyen pondéré (fractions massiques)
w = COMPOSITION
α_mean = (
    w.C3S .* α_C3S .+ w.C2S .* α_C2S .+
        w.C3A .* α_C3A .+ w.C4AF .* α_C4AF
) ./
    (w.C3S + w.C2S + w.C3A + w.C4AF)

# Résumé
println()
println("╔══════════════════════════════════════════════════╗")
println("║   Cinétique clinker — KineticsProblem            ║")
println("║   CEM I, e/c = $WC, T₀ = 20 °C, 7 jours         ║")
println("╠══════════════════════════════════════════════════╣")
@printf "║  ΔT max      = %6.2f °C\n"  maximum(T_°C_vec) - 20.0
@printf "║  α(C3S)      = %6.4f\n"    α_C3S[end]
@printf "║  α(C2S)      = %6.4f\n"    α_C2S[end]
@printf "║  α(C3A)      = %6.4f\n"    α_C3A[end]
@printf "║  α(C4AF)     = %6.4f\n"    α_C4AF[end]
@printf "║  ᾱ moyen      = %6.4f\n"    α_mean[end]
@printf "║  Q total     = %6.1f kJ/kg\n" Q_kJ_vec[end]
println("╚══════════════════════════════════════════════════╝")

# ── 8. Tracés ─────────────────────────────────────────────────────────────────

using Plots
gr()

p1 = plot(
    t_T ./ 3600, T_°C_vec;
    xlabel = "Temps [h]", ylabel = "T [°C]",
    title = "Température (calorimètre semi-adiabatique)",
    label = "T(t)", lw = 2, color = :red,
)
hline!(p1, [20.0]; linestyle = :dash, color = :gray, label = "T₀ = T_env")

p2 = plot(
    t_h, [α_C3S α_C2S α_C3A α_C4AF α_mean];
    xlabel = "Temps [h]", ylabel = "Degré d'hydratation α",
    title = "Hydratation des phases clinker",
    label = ["C₃S" "C₂S" "C₃A" "C₄AF" "ᾱ moyen"],
    lw = 2, ls = [:solid :dash :dot :dashdot :solid],
)
hline!(p2, [α_max]; linestyle = :dash, color = :black, label = "α_max (Powers)")

p3 = plot(
    t_Q ./ 3600, Q_kJ_vec;
    xlabel = "Temps [h]", ylabel = "Q [kJ/kg ciment]",
    title = "Chaleur cumulée",
    label = "Q(t)", lw = 2, color = :purple,
)

display(
    plot(
        p1, p2, p3; layout = (1, 3), size = (1400, 420),
        plot_title = "ChemistryLab — KineticsProblem cement — CEM I e/c=$WC",
    )
)
