# =============================================================================
# lavergne2018_cement_calorimetry.jl
#
# Simulation du dégagement de chaleur d'une pâte de CEM I en calorimètre
# semi-adiabatique par la méthode de Lavergne et al. (2018).
#
# Cinétique d'hydratation : modèle de Parrot & Killoh (1984)
#   via ChemistryLab.ParrotKillohRateModel (PK_PARAMS_C3S, etc.),
#   température : correction d'Arrhenius par phase,
#   eau disponible : limite maximale d'hydratation α_max (Powers 1948).
#
# Calorimètre : semi-adiabatique, perte de chaleur quadratique
#   φ(ΔT) = a·ΔT + b·ΔT²
#   (ChemistryLab.SemiAdiabaticCalorimeter avec heat_loss fonctionnel)
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
#
# Packages requis (à ajouter si absents) :
#   julia> using Pkg
#   julia> Pkg.activate("scripts")
#   julia> Pkg.add(["OrdinaryDiffEq", "Plots"])
# =============================================================================

using Pkg
Pkg.activate(@__DIR__)

using ChemistryLab
using OrdinaryDiffEq

# ── 1. Composition de la pâte de ciment ──────────────────────────────────────
#
# CEM I 52,5 type, d'après Lavergne et al. (2018)
# Fractions massiques dans le clinker (somme ≈ 0.951 ; reste = gypse, calcaire)

const PHASE_MASS_FRAC = (C3S = 0.619, C2S = 0.165, C3A = 0.080, C4AF = 0.087)
const MOLAR_MASS      = (C3S = 228.33, C2S = 172.25, C3A = 270.19, C4AF = 485.96) # g/mol

# Rapport eau/ciment
const WC = 0.40

# Moles de chaque phase par kg de ciment [mol/kg_ciment]
const N0 = map(
    (f, M) -> 1000.0 * f / M,
    PHASE_MASS_FRAC,
    MOLAR_MASS,
)  # NamedTuple (C3S=..., C2S=..., C3A=..., C4AF=...)

# Degré maximal d'hydratation (loi de Powers 1948 : α_max ≤ w/c / 0.42)
const ALPHA_MAX = min(1.0, WC / 0.42)

# ── 2. Modèles cinétiques de Parrot & Killoh depuis la bibliothèque ─────────
#
# ChemistryLab exporte PK_PARAMS_C3S / C2S / C3A / C4AF avec α_max = 1.0.
# On redéfinit les modèles avec le α_max calculé depuis le rapport e/c
# (loi de Powers 1948) pour que la cinétique s'arrête à la bonne valeur.
#
# Paramètres K₁, N₁, K₂, N₂, K₃, N₃, B, Ea — Parrot & Killoh (1984)
# Ea — Schindler & Folliard (2005) / van Breugel (1991) [J/mol]

const PK_C3S  = ParrotKillohRateModel(1.5,   3.3,  0.018,   2.5, 0.0024, 4.0, 0.5,  41_570.0; α_max = ALPHA_MAX)
const PK_C2S  = ParrotKillohRateModel(0.95,  0.5,  0.0005,  2.5, 0.0024, 4.0, 0.2,  43_670.0; α_max = ALPHA_MAX)
const PK_C3A  = ParrotKillohRateModel(0.082, 0.87, 0.00024, 2.0, 0.0024, 4.0, 0.04, 54_040.0; α_max = ALPHA_MAX) # avec sulfate
const PK_C4AF = ParrotKillohRateModel(0.165, 3.7,  0.0015,  2.5, 0.0024, 4.0, 0.5,  34_420.0; α_max = ALPHA_MAX)

# ── 3. Enthalpies d'hydratation [J/mol] ──────────────────────────────────────
#
# Valeurs en J/g d'après Taylor (1997) "Cement Chemistry", converties en J/mol.
# C3A : valeur pour l'hydratation avec gypse → C₃A + 3C̄SH₂ + 26H → C₆AŜ₃H₃₂
# (ettringite), plus exothermique que l'aluminate seul.

const ENTHALPY_J_PER_G = (C3S = -517.0, C2S = -262.0, C3A = -1144.0, C4AF = -418.0)
const ENTHALPY         = map(
    (h, M) -> h * M,
    ENTHALPY_J_PER_G,
    MOLAR_MASS,
) # J/mol

# ── 4. Calorimètre semi-adiabatique ──────────────────────────────────────────
#
# Dispositif de type Langavant (norme NF EN 196-9).
# Perte de chaleur quadratique φ(ΔT) = a·ΔT + b·ΔT² (Lavergne et al. 2018).
# Les paramètres a et b sont typiquement calibrés depuis un essai de référence ;
# les valeurs ci-dessous sont représentatives d'un vase Dewar de 2 L.
#
# Capacité thermique Cp [J/K] pour 1 kg ciment + WC kg eau + vase :
#   cp_ciment ≈ 800 J/(kg·K), cp_eau = 4186 J/(kg·K), cp_vase ≈ 900 J/(kg·K)

const T0_K    = 20.0 + 273.15   # température initiale [K]
const T_ENV_K = 20.0 + 273.15   # température ambiante [K]

const CP = 1.0 * 800.0 + WC * 4186.0 + 1.0 * 900.0   # ≈ 3449 J/K

const A_CAL = 0.30   # coefficient linéaire de perte de chaleur  [W/K]
const B_CAL = 0.003  # coefficient quadratique                   [W/K²]

cal = SemiAdiabaticCalorimeter(;
    T0        = T0_K,
    T_env     = T_ENV_K,
    Cp        = CP,
    heat_loss = ΔT -> A_CAL * ΔT + B_CAL * ΔT^2,
)

# ── 5. Définition de l'ODE ────────────────────────────────────────────────────
#
# Vecteur d'état : u = [α_C3S, α_C2S, α_C3A, α_C4AF, T]
#   dα_k/dt = r_k / n₀_k    (r_k en mol/s depuis ParrotKillohRateModel)
#   dT/dt   = (q̇ - φ(T - T_env)) / Cp
#
# Flux de chaleur instantané [W/kg_ciment] :
#   q̇ = -Σ_k r_k · ΔH_k   (> 0 car ΔH < 0 et r > 0)
#
# Note : ce script construit directement un ODEProblem avec les modèles de
# cinétique de ChemistryLab.ParrotKillohRateModel, sans passer par
# KineticsProblem, car les phases clinker ne sont pas encore dans la base de
# données thermodynamique ChemistryLab.  Le vecteur d'état est α (degré
# d'hydratation) plutôt que n (moles) pour rester cohérent avec Lavergne et al.
# La conversion entre les deux est : n_current = n₀ · (1 - α).

function cement_ode!(du, u, p, _t)
    α_C3S, α_C2S, α_C3A, α_C4AF, T = u
    pk_C3S, pk_C2S, pk_C3A, pk_C4AF, n0, ΔH, cal_local = p

    # Taux de dissolution [mol/s] via ParrotKillohRateModel
    # n_current = n₀ · (1 - α), n_initial = n₀
    r_C3S  = pk_C3S( ; T, n_current = n0.C3S  * (1 - α_C3S),  n_initial = n0.C3S)
    r_C2S  = pk_C2S( ; T, n_current = n0.C2S  * (1 - α_C2S),  n_initial = n0.C2S)
    r_C3A  = pk_C3A( ; T, n_current = n0.C3A  * (1 - α_C3A),  n_initial = n0.C3A)
    r_C4AF = pk_C4AF(; T, n_current = n0.C4AF * (1 - α_C4AF), n_initial = n0.C4AF)

    # dα/dt = r_k / n₀_k  [s⁻¹]
    du[1] = r_C3S  / n0.C3S
    du[2] = r_C2S  / n0.C2S
    du[3] = r_C3A  / n0.C3A
    du[4] = r_C4AF / n0.C4AF

    # Flux de chaleur [W/kg_ciment] : ΔH < 0 → le produit r·(-ΔH) > 0
    qdot = -(r_C3S * ΔH.C3S + r_C2S * ΔH.C2S + r_C3A * ΔH.C3A + r_C4AF * ΔH.C4AF)

    ΔT = T - cal_local.T_env
    du[5] = (qdot - cal_local.heat_loss(ΔT)) / cal_local.Cp

    return nothing
end

# ── 6. Résolution ─────────────────────────────────────────────────────────────
#
# Conditions initiales : α₀ ≈ 0 pour toutes les phases, T₀ = température initiale.
# Un petit α₀ > 0 évite l'indétermination numérique du terme de diffusion PK.
# Rodas5P est un solveur de Rosenbrock d'ordre 5 adapté aux systèmes raides.

const α₀ = 1.0e-8

u0 = [α₀, α₀, α₀, α₀, T0_K]
tspan = (0.0, 7.0 * 86400.0)   # 7 jours [s]

p_ode = (
    pk_C3S  = PK_C3S,
    pk_C2S  = PK_C2S,
    pk_C3A  = PK_C3A,
    pk_C4AF = PK_C4AF,
    n0      = N0,
    ΔH      = ENTHALPY,
    cal     = cal,
)

prob = ODEProblem(cement_ode!, u0, tspan, p_ode)

@info "Résolution en cours (7 jours, Rodas5P)..."
sol = solve(prob, Rodas5P(); reltol = 1.0e-6, abstol = 1.0e-9, saveat = 900.0)
@info "Résolution terminée : $(length(sol.t)) points de sauvegarde."

# ── 7. Post-traitement ────────────────────────────────────────────────────────

t_h = sol.t ./ 3600.0   # secondes → heures

T_K_vec  = [u[5] for u in sol.u]
T_°C_vec = T_K_vec .- 273.15

α_C3S_vec  = [u[1] for u in sol.u]
α_C2S_vec  = [u[2] for u in sol.u]
α_C3A_vec  = [u[3] for u in sol.u]
α_C4AF_vec = [u[4] for u in sol.u]

# ── Flux de chaleur q̇(t) récupéré via le bilan du calorimètre ──────────────
# q̇ = Cp·dT/dt + φ(T - T_env)   [W/kg_ciment]
qdot_W = zeros(length(sol.t))
for i in 2:lastindex(sol.t)
    dt   = sol.t[i] - sol.t[i - 1]
    dTdt = (T_K_vec[i] - T_K_vec[i - 1]) / dt
    ΔT   = T_K_vec[i] - T_ENV_K
    qdot_W[i] = CP * dTdt + cal.heat_loss(ΔT)
end

# ── Chaleur cumulée Q(t) = ∫₀ᵗ q̇ dτ  [J/kg_ciment → kJ/kg_ciment] ─────────
Δt_vec = vcat(0.0, diff(sol.t))
Q_kJ   = cumsum(qdot_W .* Δt_vec) ./ 1000.0

# ── Degré d'hydratation moyen pondéré par la fraction massique ───────────────
α_mean = (PHASE_MASS_FRAC.C3S .* α_C3S_vec  .+
          PHASE_MASS_FRAC.C2S .* α_C2S_vec  .+
          PHASE_MASS_FRAC.C3A .* α_C3A_vec  .+
          PHASE_MASS_FRAC.C4AF .* α_C4AF_vec)

# ── Résumé ────────────────────────────────────────────────────────────────────
println()
println("╔══════════════════════════════════════════════╗")
println("║    Résultats calorimétrie semi-adiabatique   ║")
println("║    CEM I — Lavergne et al. (2018)            ║")
println("╠══════════════════════════════════════════════╣")
println("║  Temps       = 7 jours")
println("║  T final     = $(round(T_°C_vec[end], digits=2)) °C")
println("║  ΔT max      = $(round(maximum(T_°C_vec) - (T0_K - 273.15), digits=2)) °C")
println("║  α(C3S)      = $(round(α_C3S_vec[end],  sigdigits=3))")
println("║  α(C2S)      = $(round(α_C2S_vec[end],  sigdigits=3))")
println("║  α(C3A)      = $(round(α_C3A_vec[end],  sigdigits=3))")
println("║  α(C4AF)     = $(round(α_C4AF_vec[end], sigdigits=3))")
println("║  α moyen     = $(round(α_mean[end], sigdigits=3))")
println("║  Q total     = $(round(Q_kJ[end], digits=1)) kJ/kg_ciment")
println("╚══════════════════════════════════════════════╝")

# ── 8. Tracés (nécessite Plots) ───────────────────────────────────────────────
# Décommenter si Plots est disponible dans l'environnement :
#
using Plots
gr()

p1 = plot(t_h, T_°C_vec;
    xlabel = "Temps [h]", ylabel = "Température [°C]",
    title  = "Profil de température — calorimètre semi-adiabatique",
    label  = "T(t)", lw = 2, color = :red)
hline!(p1, [T0_K - 273.15]; linestyle = :dash, color = :gray, label = "T₀ = T_env")

p2 = plot(t_h, [α_C3S_vec α_C2S_vec α_C3A_vec α_C4AF_vec α_mean];
    xlabel = "Temps [h]", ylabel = "Degré d'hydratation α",
    title  = "Hydratation des phases clinker (modèle Parrot-Killoh)",
    label  = ["C₃S" "C₂S" "C₃A" "C₄AF" "ᾱ moyen"],
    lw = 2, ls = [:solid :dash :dot :dashdot :solid])
hline!(p2, [ALPHA_MAX]; linestyle = :dash, color = :black, label = "α_max")

p3 = plot(t_h, qdot_W ./ 1000;
    xlabel = "Temps [h]", ylabel = "Flux de chaleur [kW/kg_ciment]",
    title  = "Flux thermique instantané",
    label  = "q̇(t)", lw = 2, color = :orange)

p4 = plot(t_h, Q_kJ;
    xlabel = "Temps [h]", ylabel = "Q [kJ/kg_ciment]",
    title  = "Chaleur cumulée",
    label  = "Q(t)", lw = 2, color = :purple)

display(plot(p1, p2, p3, p4; layout = (2, 2), size = (1200, 800),
    plot_title = "Lavergne et al. (2018) — CEM I, E/C = $WC"))
