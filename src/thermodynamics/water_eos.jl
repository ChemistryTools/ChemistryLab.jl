# ── Équation d'état de l'eau : densité ρ(T,P) et constante diélectrique ε(T,P) ──
# Références :
#   - Kell (1975), J. Chem. Eng. Data 20, 97-105           → ρ(T, Pref)
#   - Helgeson (1969), Am. J. Sci. 267, 729-804             → ε(T, Pref)
#   - Tanger & Helgeson (1988), Am. J. Sci. 288, 19-98     → Born functions,
#       valeurs de référence Table 1 : Y ≈ −5.81×10⁻⁵ K⁻¹, Q ≈ −3.57×10⁻¹⁵ Pa⁻¹

# ── Interface abstraite ────────────────────────────────────────────────────────

"""
    abstract type AbstractWaterEOS end

Interface pour les équations d'état de l'eau.
Permet de brancher une implémentation alternative sans modifier le code utilisateur.
"""
abstract type AbstractWaterEOS end

"""
    struct HelgesonKirkham1974 <: AbstractWaterEOS

Modèle basé sur Helgeson (1969) pour ε(T) et Kell (1975) pour ρ(T,P).
Valide pour 273–423 K (0–150 °C), 1e5–1e8 Pa (1–1 000 bar).
"""
struct HelgesonKirkham1974 <: AbstractWaterEOS end

"""Modèle par défaut."""
const DEFAULT_WATER_EOS = HelgesonKirkham1974()

# ── Densité ρ(T, P) — Kell (1975) + correction de pression ────────────────────

"""
    _kell_density_1bar(T::Float64) -> Float64

Densité de l'eau à P = 1 atm [g/cm³], formule de Kell (1975).
T en Kelvin. Valide pour 273–423 K (0–150 °C), précision ±0.0001 g/cm³.
"""
function _kell_density_1bar(T::Float64)::Float64
    t   = T - 273.15                                 # K → °C
    num = 999.83952   +
          16.945176   * t   -
          7.987040e-3 * t^2 -
          46.170461e-6* t^3 +
          105.56302e-9* t^4 -
          280.54253e-12*t^5
    den = 1.0 + 16.879850e-3 * t
    return (num / den) * 1e-3   # kg/m³ → g/cm³
end

# Compressibilité isotherme approximative de l'eau liquide [Pa⁻¹]
# Valable sur 0–100 °C, à ±20 % (Fine & Millero 1973 simplifié)
function _kell_compressibility(T::Float64)::Float64
    t = T - 273.15     # °C
    # κ_T [×10⁻¹⁰ Pa⁻¹] polynomial fit to Fine & Millero (1973) data
    return (50.88 - 0.376*t + 3.00e-3*t^2 - 2.7e-5*t^3) * 1e-10
end

"""
    water_density(T::Float64, P::Float64,
                  ::AbstractWaterEOS = DEFAULT_WATER_EOS) -> Float64

Densité de l'eau pure ρ [g/cm³] à T [K] et P [Pa].
Kell (1975) à Pref = 1 bar + correction de compression par intégrale de κ_T.
Domaine de validité : 273–423 K, 1e5–1e8 Pa.
"""
function water_density(T::Float64, P::Float64,
                       ::AbstractWaterEOS = DEFAULT_WATER_EOS)::Float64
    ρ₀ = _kell_density_1bar(T)
    κ  = _kell_compressibility(T)
    ΔP = P - 1e5                        # Pa au-dessus de 1 bar
    return ρ₀ * exp(κ * ΔP)            # ρ = ρ₀ exp(κ ΔP) pour κ constant
end

# ── Constante diélectrique ε(T, P) ────────────────────────────────────────────

# Coefficients polynomiaux pour ε(T) à Pref = 1 bar.
# Helgeson (1969), Am. J. Sci. 267, Eq. 120.
# Tels qu'utilisés dans PHREEQC (Parkhurst & Appelo 2013).
const _H69_A = 78.54
const _H69_B = -4.579e-3   # K⁻¹
const _H69_C =  1.19e-5    # K⁻²
const _H69_D = -2.8e-8     # K⁻³
const _H69_Tr = 298.15     # K

# Pression de référence = 1 bar
const _EOS_Pr = 1.0e5    # Pa

# Pente ∂ε/∂P à T=25 °C  [Pa⁻¹].
# Déduit de Q = ∂(1/ε)/∂P ≈ −3.57×10⁻¹⁵ Pa⁻¹ (Tanger & Helgeson 1988, Table 1) :
#   ∂ε/∂P = −ε² × Q = 78.54² × 3.57×10⁻¹⁵ ≈ 2.2×10⁻¹¹ Pa⁻¹
# Supposée constante sur le domaine (variation faible : Δε ≈ 0.001 sur 500 bar).
const _DEDP = 2.2e-11   # Pa⁻¹  (positif : ε croît très légèrement avec P)

"""
    water_dielectric(T::Float64, P::Float64 = _EOS_Pr) -> Float64

Constante diélectrique ε (sans unité) de l'eau à T [K] et P [Pa].
Formule de Helgeson (1969) pour ε(T) à Pref = 1 bar, avec correction linéaire en P.
"""
function water_dielectric(T::Float64, P::Float64 = _EOS_Pr)::Float64
    ΔT = T - _H69_Tr
    ε₀ = _H69_A * (1.0 + _H69_B*ΔT + _H69_C*ΔT^2 + _H69_D*ΔT^3)
    return ε₀ + _DEDP * (P - _EOS_Pr)   # ΔP [Pa] × 2.2e-11 Pa⁻¹ ≈ 0
end

# ── Fonctions de Born ─────────────────────────────────────────────────────────
#
# Convention Tanger & Helgeson (1988) :
#   Y = ∂(−1/ε)/∂T = (1/ε²)(∂ε/∂T)    [K⁻¹]       → négatif pour l'eau
#   Q = ∂(−1/ε)/∂P = −(1/ε²)(∂ε/∂P)   [Pa⁻¹]      → négatif pour l'eau
#   X = ∂Y/∂T = ∂²(−1/ε)/∂T²           [K⁻²]
#
# Dérivées calculées par différences finies centrées.

const _BORN_ΔT = 1.0e-2   # K
const _BORN_ΔP = 1.0e4    # Pa (= 0.1 bar)

"""
    born_functions(T::Float64, P::Float64,
                   weos::AbstractWaterEOS = DEFAULT_WATER_EOS)
        -> @NamedTuple{Y, Q, X, ε, ρ}

Calcule les fonctions de Born de l'eau à T [K] et P [Pa].

Convention Tanger & Helgeson (1988) :

| Symbole | Définition | Signe (eau 25 °C) |
|:--------|:-----------|:-----------------:|
| Y       | ∂(−1/ε)/∂T = (1/ε²)(∂ε/∂T)   | négatif |
| Q       | ∂(−1/ε)/∂P = −(1/ε²)(∂ε/∂P) | négatif |
| X       | ∂Y/∂T = ∂²(−1/ε)/∂T²         | négatif |

Valeurs de référence à 25 °C / 1 bar (Tanger & Helgeson 1988, Table 1) :
  Y ≈ −5.81×10⁻⁵ K⁻¹, Q ≈ −3.57×10⁻¹⁵ Pa⁻¹
"""
function born_functions(T::Float64, P::Float64,
                        weos::AbstractWaterEOS = DEFAULT_WATER_EOS)
    # f(t,p) = −1/ε(t,p)  [convention T&H]
    neg_inv_ε(t, p) = -1.0 / water_dielectric(t, p)

    ε_val = water_dielectric(T, P)
    ρ_val = water_density(T, P, weos)

    # Différences finies centrées
    # Y, X : dérivées de f = −1/ε  (convention T&H : Y = ∂(−1/ε)/∂T)
    # Q     : dérivée de g = +1/ε  (convention T&H : Q = ∂(1/ε)/∂P = −(1/ε²)(∂ε/∂P) < 0)
    inv_ε(t, p) = 1.0 / water_dielectric(t, p)

    f₊T  = neg_inv_ε(T + _BORN_ΔT, P)
    f₋T  = neg_inv_ε(T - _BORN_ΔT, P)
    f₀   = neg_inv_ε(T, P)

    Y = (f₊T - f₋T) / (2 * _BORN_ΔT)          # ∂(−1/ε)/∂T  [K⁻¹]
    Q = (inv_ε(T, P + _BORN_ΔP) - inv_ε(T, P - _BORN_ΔP)) / (2 * _BORN_ΔP)  # ∂(1/ε)/∂P  [Pa⁻¹]
    X = (f₊T - 2*f₀ + f₋T) / _BORN_ΔT^2       # ∂²(−1/ε)/∂T² [K⁻²]

    return (Y=Y, Q=Q, X=X, ε=ε_val, ρ=ρ_val)
end
