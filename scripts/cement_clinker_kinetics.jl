# =============================================================================
# cement_clinker_kinetics.jl
#
# Simulation cinétique de l'hydratation du clinker OPC via ChemistryLab :
#   - ChemicalSystem construit depuis la base de données CEMDATA18
#   - KineticsProblem avec ParrotKillohRateModel (4 phases clinker)
#   - Suivi de la chaleur via SemiAdiabaticCalorimeter
#
# Ce script montre l'utilisation de la chaîne complète :
#   1. ChemicalSystem  → 2. ChemicalState  → 3. KineticReaction
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
using Printf

# ── 1. ChemicalSystem depuis la base CEMDATA18 ────────────────────────────────
#
# On sélectionne :
#   • les 4 phases clinker (espèces cinétiques)
#   • les principaux produits d'hydratation (équilibre instantané)
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

# Amorçage pH neutre pour H⁺ / OH⁻ (si présents dans le système)
try
    V0 = volume(state0)
    set_quantity!(state0, "H+", 1.0e-7u"mol/L" * V0.liquid)
    set_quantity!(state0, "OH-", 1.0e-7u"mol/L" * V0.liquid)
catch
    # espèces absentes du système sélectionné — OK
end

# ── 3. Modèles cinétiques de Parrot & Killoh ──────────────────────────────────
#
# Degré maximal d'hydratation selon Powers (1948) : α_max ≤ w/c / 0.42
# On crée des modèles avec le α_max calculé plutôt que la valeur par défaut 1.0.

const α_max = min(1.0, WC / 0.42)

function _pk_with_amax(template::ParrotKillohRateModel, α::Real)
    return ParrotKillohRateModel(
        template.K₁, template.N₁, template.K₂, template.N₂,
        template.K₃, template.N₃, template.B, template.Ea;
        T_ref = template.T_ref, α_max = α,
    )
end

pk_C3S = _pk_with_amax(PK_PARAMS_C3S, α_max)
pk_C2S = _pk_with_amax(PK_PARAMS_C2S, α_max)
pk_C3A = _pk_with_amax(PK_PARAMS_C3A, α_max)
pk_C4AF = _pk_with_amax(PK_PARAMS_C4AF, α_max)

# ── 4. Réactions cinétiques ───────────────────────────────────────────────────
#
# Le constructeur de commodité KineticReaction(cs, nom, modèle, surface)
# recherche automatiquement l'indice de l'espèce dans cs et construit
# la stœchiométrie par défaut (-1 pour le minéral, 0 pour le reste).
# FixedSurfaceArea(1.0) est un dummy car ParrotKillohRateModel n'utilise
# pas A_surface.
#
# heat_per_mol : chaleurs d'hydratation standard [J/mol] (Taylor 1997 / Lothenbach 2006)
#   C3S : 502 J/g × 228.32 g/mol ≈ 114 617 J/mol
#   C2S : 260 J/g × 172.24 g/mol ≈  44 782 J/mol
#   C3A : 866 J/g × 270.19 g/mol ≈ 233 985 J/mol
#   C4AF: 419 J/g × 485.96 g/mol ≈ 203 617 J/mol

kr_C3S = KineticReaction(cs, "C3S", pk_C3S, FixedSurfaceArea(1.0); heat_per_mol = 114_617.0)
kr_C2S = KineticReaction(cs, "C2S", pk_C2S, FixedSurfaceArea(1.0); heat_per_mol = 44_782.0)
kr_C3A = KineticReaction(cs, "C3A", pk_C3A, FixedSurfaceArea(1.0); heat_per_mol = 233_985.0)
kr_C4AF = KineticReaction(cs, "C4AF", pk_C4AF, FixedSurfaceArea(1.0); heat_per_mol = 203_617.0)

@info "Indice C3S dans cs : $(kr_C3S.idx_mineral)"

# ── 5. Problème cinétique + calorimètre ──────────────────────────────────────
#
# tspan : 7 jours en secondes
# equilibrium_solver = nothing → pas de re-spéciation aqueuse (suffisant pour PK)

const TSPAN = (0.0, 7.0 * 86400.0)

kp = KineticsProblem(
    cs,
    [kr_C3S, kr_C2S, kr_C3A, kr_C4AF],
    state0,
    TSPAN;
    equilibrium_solver = nothing,
)

# Paramètres calorimètre semi-adiabatique (Lavergne et al. 2018)
#   Cp [J/K] : 1 kg ciment + WC kg eau + vase Dewar
#   pertes quadratiques : φ(ΔT) = a·ΔT + b·ΔT²
cal = SemiAdiabaticCalorimeter(;
    T0 = 293.15,     # température initiale [K]  (20 °C)
    T_env = 293.15,     # température ambiante [K]
    Cp = 1.0 * 800.0 + WC * 4186.0 + 1.0 * 900.0,   # ≈ 3449 J/K
    heat_loss = ΔT -> 0.3 * ΔT + 0.003 * ΔT^2,
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

# Degrés d'hydratation depuis les moles cinétiques
n0_kin = sol.prob.p.n_initial     # moles initiales [mol] de chaque phase cinétique
n_kin = [[u[i] for u in sol.u] for i in eachindex(n0_kin)]

# Indices des phases dans le vecteur d'état ODE (ordre = [kr_C3S, kr_C2S, kr_C3A, kr_C4AF])
α_C3S = 1.0 .- n_kin[1] ./ n0_kin[1]
α_C2S = 1.0 .- n_kin[2] ./ n0_kin[2]
α_C3A = 1.0 .- n_kin[3] ./ n0_kin[3]
α_C4AF = 1.0 .- n_kin[4] ./ n0_kin[4]

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
    label = "T(t)", lw = 2, color = :red
)
hline!(p1, [20.0]; linestyle = :dash, color = :gray, label = "T₀ = T_env")

p2 = plot(
    t_h, [α_C3S α_C2S α_C3A α_C4AF α_mean];
    xlabel = "Temps [h]", ylabel = "Degré d'hydratation α",
    title = "Hydratation des phases clinker",
    label = ["C₃S" "C₂S" "C₃A" "C₄AF" "ᾱ moyen"],
    lw = 2, ls = [:solid :dash :dot :dashdot :solid]
)
hline!(p2, [α_max]; linestyle = :dash, color = :black, label = "α_max (Powers)")

p3 = plot(
    t_Q ./ 3600, Q_kJ_vec;
    xlabel = "Temps [h]", ylabel = "Q [kJ/kg ciment]",
    title = "Chaleur cumulée",
    label = "Q(t)", lw = 2, color = :purple
)

display(
    plot(
        p1, p2, p3; layout = (1, 3), size = (1400, 420),
        plot_title = "ChemistryLab — KineticsProblem cement — CEM I e/c=$WC"
    )
)
