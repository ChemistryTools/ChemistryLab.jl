# =============================================================================
# blended_cement_kinetics.jl
#
# Simulation cinétique de l'hydratation d'un ciment avec additions minérales :
#   - 63% clinker OPC (CEM I 52.5 R, 4 phases Parrot–Killoh)
#   - 30% laitier de haut-fourneau (GGBS, lent)
#   - 7% métakaolin (MK, plus réactif que le laitier)
#
# CEMDATA18 ne contient pas GGBS ni MK comme réactants : on crée des Species
# personnalisées avec un ΔₐG⁰ factice (le modèle Parrot–Killoh ignore Ω).
# La mole de référence de chaque addition est définie par le champ :M explicite
# de la Species (masse par « unité formulaire représentative »), qui doit être
# cohérent avec heat_per_mol pour que la calorimétrie soit correcte :
#
#     heat_per_mol [J/mol] = chaleur_spec [J/g] × M [g/mol]
#
# Usage :
#   julia --project scripts/blended_cement_kinetics.jl
#   ou depuis le REPL : include("scripts/blended_cement_kinetics.jl")
# =============================================================================

using Pkg
Pkg.activate(@__DIR__)

using ChemistryLab
using OrdinaryDiffEq
using DynamicQuantities
using Printf

# ── 1. ChemicalSystem de base depuis CEMDATA18 ───────────────────────────────
#
# On sélectionne les phases clinker, les principaux produits d'hydratation
# (portlandite, C-S-H, ettringite, AFm, hydrotalcite, strätlingite…) et l'eau.
# GGBS et MK seront ajoutés ensuite en tant qu'espèces personnalisées.

const DATA_FILE = joinpath(pkgdir(ChemistryLab), "data", "cemdata18-thermofun.json")

substances = build_species(DATA_FILE)

input_species = split(
    "C3S C2S C3A C4AF " *
        "Portlandite Jennite ettringite monosulphate12 C3AH6 C3FH6 " *
        "straetlingite hydrotalcite " *
        "H2O@",
)

species_base = speciation(substances, input_species; aggregate_state = [AS_AQUEOUS])
cs_base = ChemicalSystem(species_base, CEMDATA_PRIMARIES)

@info "Système de base : $(length(cs_base.species)) espèces"

# ── 2. Species personnalisées pour GGBS et MK ─────────────────────────────────
#
# GGBS (laitier granulé de haut-fourneau) :
#   Composition typique CEM I-like : ~42% CaO, 35% SiO₂, 12% Al₂O₃, 8% MgO.
#   Formule représentative (charge neutre) : CaAl₂Si₂O₈ (analogie anorthite).
#   Masse molaire de l'unité formulaire écrasée à 95 g/mol pour que :
#       heat_per_mol = 380 J/g × 95 g/mol ≈ 36 100 J/mol  (Gruyaert 2010)
#
# MK (métakaolin) :
#   Formule exacte : Al₂Si₂O₇ (kaolinite déshydroxylée).
#   Masse molaire auto-calculée : M = 222 g/mol.
#   heat_per_mol = 250 J/g × 222 g/mol ≈ 55 500 J/mol  (Lothenbach 2011)
#
# ΔₐG⁰ factice : très négatif → Ω ≈ 0 (dissolution toujours favorisée).
# Parrot–Killoh ignore Ω ; la valeur n'influence pas les taux cinétiques.

const _dummy_G = NumericFunc((T, P) -> -1_200_000.0, (:T, :P), u"J/mol")

# GGBS : formule CaAl₂Si₂O₈, :M écrasé à 95 g/mol
sp_ggbs = Species(
    "CaAl2Si2O8";
    symbol = "GGBS",
    name = "GGBS",
    aggregate_state = AS_CRYSTAL,
    properties = Dict{Symbol, Any}(
        :M => 0.095u"kg/mol",
        :ΔₐG⁰ => _dummy_G,
    ),
)

# MK : formule Al₂Si₂O₇ (métakaolin), M = 222 g/mol auto-calculé
sp_mk = Species(
    "Al2Si2O7";
    symbol = "MK",
    name = "MK",
    aggregate_state = AS_CRYSTAL,
    properties = Dict{Symbol, Any}(
        :ΔₐG⁰ => _dummy_G,
    ),
)

# ── 3. ChemicalSystem étendu ──────────────────────────────────────────────────

all_species = vcat(cs_base.species, sp_ggbs, sp_mk)
cs = ChemicalSystem(all_species, CEMDATA_PRIMARIES)

@info "Système étendu : $(length(cs.species)) espèces (dont GGBS et MK)"

# ── 4. État initial ───────────────────────────────────────────────────────────
#
# Ciment ternaire (1 kg) : 63% clinker OPC, 30% laitier GGBS, 7% métakaolin
# Ratio eau/liant w/b = 0.40
#
# Phases clinker (fractions massiques dans le clinker CEM I 52.5 R) :
#   C₃S : 61.9%, C₂S : 16.5%, C₃A : 8.0%, C₄AF : 8.7%   (Lavergne 2018)
# Fraction clinker dans le ciment ternaire : 0.63
#
# heat_per_mol [J/mol] (par unité formulaire définie par :M) :
#   C₃S  : 502 J/g × 228 g/mol ≈ 114 617  (Taylor 1997)
#   C₂S  : 260 J/g × 172 g/mol ≈  44 782  (Taylor 1997)
#   C₃A  : 866 J/g × 270 g/mol ≈ 233 985  (Taylor 1997)
#   C₄AF : 419 J/g × 486 g/mol ≈ 203 617  (Taylor 1997)
#   GGBS : 380 J/g ×  95 g/mol ≈  36 100  (Gruyaert 2010)
#   MK   : 250 J/g × 222 g/mol ≈  55 500  (Lothenbach 2011)

const WB = 0.4      # rapport eau/liant

const CLINKER_FRAC = 0.63      # fraction massique du clinker dans le ciment
const GGBS_FRAC = 0.3      # fraction massique du laitier
const MK_FRAC = 0.07      # fraction massique du métakaolin

const CLINKER_COMP = (C3S = 0.619, C2S = 0.165, C3A = 0.08, C4AF = 0.087)

const COMPOSITION = (
    C3S = CLINKER_COMP.C3S * CLINKER_FRAC,
    C2S = CLINKER_COMP.C2S * CLINKER_FRAC,
    C3A = CLINKER_COMP.C3A * CLINKER_FRAC,
    C4AF = CLINKER_COMP.C4AF * CLINKER_FRAC,
    GGBS = GGBS_FRAC,
    MK = MK_FRAC,
)

state0 = ChemicalState(cs)
for (name, frac) in pairs(COMPOSITION)
    set_quantity!(state0, string(name), frac * u"kg")
end
set_quantity!(state0, "H2O@", WB * u"kg")

# ── 5. Modèles cinétiques ─────────────────────────────────────────────────────
#
# α_max selon Powers (1948) : hydratation limitée par l'eau disponible
const α_max = min(1.0, WB / 0.42)

function _pk_with_amax(template::ParrotKillohRateModel, α::Real)
    return ParrotKillohRateModel(
        template.K₁, template.N₁, template.K₂, template.N₂,
        template.K₃, template.N₃, template.B, template.Ea;
        T_ref = template.T_ref, α_max = α,
    )
end

# Clinker (Parrot & Killoh 1984, corrections Schindler & Folliard 2005)
pk_C3S = _pk_with_amax(PK_PARAMS_C3S, α_max)
pk_C2S = _pk_with_amax(PK_PARAMS_C2S, α_max)
pk_C3A = _pk_with_amax(PK_PARAMS_C3A, α_max)
pk_C4AF = _pk_with_amax(PK_PARAMS_C4AF, α_max)

# GGBS (laitier) : paramètres PK adaptés de la littérature
#   K₁ faible → cinétique initiale lente (vitreux, peu de nucléation)
#   Ea élevée → forte sensibilité à la température
#   α_max = 0.90 (vitreux, hydratation partielle)
#   Références : Richardson & Groves (1992), Chen & Brouwers (2007)
pk_ggbs = ParrotKillohRateModel(
    0.15, 2.0, 0.003, 2.0, 0.0015, 3.5, 0.2, 46_000.0;
    T_ref = 293.15, α_max = 0.9,
)

# MK (métakaolin) : paramètres PK adaptés de la littérature
#   K₁ élevé → réactivité initiale forte (surface spécifique élevée)
#   α_max = 0.95 (quasi-totale dans les conditions normales)
#   Références : Lothenbach et al. (2011), Deschner et al. (2012)
pk_mk = ParrotKillohRateModel(
    0.7, 1.5, 0.008, 2.0, 0.002, 3.5, 0.3, 48_000.0;
    T_ref = 293.15, α_max = 0.95,
)

# ── 6. Réactions cinétiques ───────────────────────────────────────────────────
#
# Le constructeur de commodité KineticReaction(cs, nom, modèle, surface)
# résout l'espèce par son champ :symbol dans cs.
# FixedSurfaceArea(1.0) est un dummy : Parrot–Killoh n'utilise pas A_surface.

kr_C3S = KineticReaction(
    cs, "C3S", pk_C3S, FixedSurfaceArea(1.0);
    heat_per_mol = 114_617.0,
)
kr_C2S = KineticReaction(
    cs, "C2S", pk_C2S, FixedSurfaceArea(1.0);
    heat_per_mol = 44_782.0,
)
kr_C3A = KineticReaction(
    cs, "C3A", pk_C3A, FixedSurfaceArea(1.0);
    heat_per_mol = 233_985.0,
)
kr_C4AF = KineticReaction(
    cs, "C4AF", pk_C4AF, FixedSurfaceArea(1.0);
    heat_per_mol = 203_617.0,
)
kr_ggbs = KineticReaction(
    cs, "GGBS", pk_ggbs, FixedSurfaceArea(1.0);
    heat_per_mol = 36_100.0,    # 380 J/g × 95 g/mol (Gruyaert 2010)
)
kr_mk = KineticReaction(
    cs, "MK", pk_mk, FixedSurfaceArea(1.0);
    heat_per_mol = 55_500.0,    # 250 J/g × 222 g/mol (Lothenbach 2011)
)

@info "Indices : C3S=$(kr_C3S.idx_mineral)  GGBS=$(kr_ggbs.idx_mineral)  MK=$(kr_mk.idx_mineral)"

# ── 7. Problème cinétique ─────────────────────────────────────────────────────
#
# tspan : 28 jours → montre la réactivité tardive du laitier
# equilibrium_solver = nothing : pas de re-spéciation aqueuse (suffisant pour PK)

const TSPAN = (0.0, 28.0 * 86400.0)

kp = KineticsProblem(
    cs,
    [kr_C3S, kr_C2S, kr_C3A, kr_C4AF, kr_ggbs, kr_mk],
    state0,
    TSPAN;
    equilibrium_solver = nothing,
)

# ── 8. Calorimètre semi-adiabatique ──────────────────────────────────────────
#
# Cp [J/K] : ciment ternaire + eau + vase Dewar
#   cement : 1 kg × 800 J/(kg·K)
#   eau    : WB kg × 4186 J/(kg·K)
#   Dewar  : 1 kg × 900 J/(kg·K)
# Pertes quadratiques (Lavergne et al. 2018)

cal = SemiAdiabaticCalorimeter(;
    T0 = 293.15,
    T_env = 293.15,
    Cp = 1.0 * 800.0 + WB * 4186.0 + 1.0 * 900.0,
    heat_loss = ΔT -> 0.3 * ΔT + 0.003 * ΔT^2,
)

# ── 9. Intégration ────────────────────────────────────────────────────────────

@info "Intégration en cours (28 jours, Rodas5P) …"
ks = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-9)
sol = integrate(kp, ks; calorimeter = cal)
@info "Terminé : $(length(sol.t)) pas acceptés."

# ── 10. Post-traitement ────────────────────────────────────────────────────────

# Profil température et chaleur depuis le calorimètre
t_T, T_K_vec = temperature_profile(sol, cal)
t_Q, Q_J_vec = cumulative_heat(sol, cal)
T_°C_vec = T_K_vec .- 273.15
Q_kJ_vec = Q_J_vec ./ 1000.0

# Degrés d'hydratation par phase (ordre du vecteur d'état ODE)
# ODE state : [C3S, C2S, C3A, C4AF, GGBS, MK]
n0_kin = sol.prob.p.n_initial        # moles initiales [mol]
n_kin = [[u[i] for u in sol.u] for i in eachindex(n0_kin)]
t_h = sol.t ./ 3600.0             # [s] → [h]

α_C3S = 1.0 .- n_kin[1] ./ n0_kin[1]
α_C2S = 1.0 .- n_kin[2] ./ n0_kin[2]
α_C3A = 1.0 .- n_kin[3] ./ n0_kin[3]
α_C4AF = 1.0 .- n_kin[4] ./ n0_kin[4]
α_GGBS = 1.0 .- n_kin[5] ./ n0_kin[5]
α_MK = 1.0 .- n_kin[6] ./ n0_kin[6]

# Degré d'hydratation moyen pondéré par la fraction massique de chaque phase
w = COMPOSITION
α_mean = (
    w.C3S .* α_C3S .+ w.C2S .* α_C2S .+
        w.C3A .* α_C3A .+ w.C4AF .* α_C4AF .+
        w.GGBS .* α_GGBS .+ w.MK .* α_MK
) ./ (w.C3S + w.C2S + w.C3A + w.C4AF + w.GGBS + w.MK)

# ── Résumé ────────────────────────────────────────────────────────────────────

println()
println("╔═══════════════════════════════════════════════════════════╗")
println("║   Cinétique ciment ternaire — ChemistryLab                ║")
@printf "║   63%% CEM I + 30%% GGBS + 7%% MK  │  w/b = %.2f  │  28 j   ║\n" WB
println("╠═══════════════════════════════════════════════════════════╣")
@printf "║  ΔT max      = %6.2f °C                                  ║\n" maximum(T_°C_vec) - 20.0
println("╠═══════════════════════════════════════════════════════════╣")
@printf "║  α(C₃S)  7j  = %6.4f   28j = %6.4f                    ║\n" α_C3S[findlast(t -> t <= 7 * 86400, sol.t)] α_C3S[end]
@printf "║  α(C₂S)  7j  = %6.4f   28j = %6.4f                    ║\n" α_C2S[findlast(t -> t <= 7 * 86400, sol.t)] α_C2S[end]
@printf "║  α(GGBS) 7j  = %6.4f   28j = %6.4f                    ║\n" α_GGBS[findlast(t -> t <= 7 * 86400, sol.t)] α_GGBS[end]
@printf "║  α(MK)   7j  = %6.4f   28j = %6.4f                    ║\n" α_MK[findlast(t -> t <= 7 * 86400, sol.t)] α_MK[end]
@printf "║  ᾱ moyen  7j  = %6.4f   28j = %6.4f                    ║\n" α_mean[findlast(t -> t <= 7 * 86400, sol.t)] α_mean[end]
println("╠═══════════════════════════════════════════════════════════╣")
@printf "║  Q  7j = %6.1f kJ/kg   Q 28j = %6.1f kJ/kg            ║\n" Q_kJ_vec[findlast(t -> t <= 7 * 86400, sol.t)] Q_kJ_vec[end]
println("╚═══════════════════════════════════════════════════════════╝")

# ── 11. Tracés ────────────────────────────────────────────────────────────────

using Plots
gr()

p1 = plot(
    t_T ./ 3600, T_°C_vec;
    xlabel = "Temps [h]", ylabel = "T [°C]",
    title = "Température",
    label = "T(t)", lw = 2, color = :red,
)
hline!(p1, [20.0]; linestyle = :dash, color = :gray, label = "T₀")

p2 = plot(
    t_h, [α_C3S α_C2S α_C3A α_C4AF α_GGBS α_MK α_mean];
    xlabel = "Temps [h]", ylabel = "Degré d'hydratation α",
    title = "Hydratation des phases",
    label = ["C₃S" "C₂S" "C₃A" "C₄AF" "GGBS" "MK" "ᾱ moyen"],
    lw = 2, ls = [:solid :dash :dot :dashdot :solid :dash :solid],
    color = [:blue :cyan :green :orange :brown :purple :black],
)
hline!(p2, [α_max]; linestyle = :dash, color = :black, label = "α_max")

p3 = plot(
    t_Q ./ 3600, Q_kJ_vec;
    xlabel = "Temps [h]", ylabel = "Q [kJ/kg]",
    title = "Chaleur cumulée",
    label = "Q(t)", lw = 2, color = :purple,
)

display(
    plot(
        p1, p2, p3;
        layout = (1, 3), size = (1500, 450),
        plot_title = "ChemistryLab — Ciment ternaire 63% OPC + 30% GGBS + 7% MK  (w/b=$WB)",
    )
)
