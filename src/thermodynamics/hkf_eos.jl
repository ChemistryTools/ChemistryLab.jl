# ── Équation d'état HKF révisée (Tanger & Helgeson 1988) ──────────────────────
# Référence principale : Tanger & Helgeson (1988), Am. J. Sci. 288, 19–98
# Propriétés standard Cp°, V°, S°, H°, G° des espèces aqueuses en fonction de T et P.
#
# Paramètres d'une espèce : a₁, a₂, a₃, a₄, c₁, c₂, ω (7 scalaires)
# Fonctions de Born Y, Q, X fournies par water_eos.jl

# ── Constantes universelles HKF ───────────────────────────────────────────────

"""Température singulière HKF [K] (Tanger & Helgeson 1988)."""
const HKF_Θ = 228.0        # K

"""Pression singulière HKF [Pa] (= 2600 bar)."""
const HKF_Ψ = 2.6e8        # Pa

"""Température de référence [K]."""
const HKF_Tr = 298.15      # K

"""Pression de référence [Pa] (= 1 bar)."""
const HKF_Pr = 1.0e5       # Pa

# Constante diélectrique et fonctions de Born au point de référence (25 °C, 1 bar)
# Tanger & Helgeson (1988), Table 1
const _HKF_εr = 78.47      # sans unité
const _HKF_Yr = -5.81e-5   # K⁻¹
const _HKF_Qr = -3.57e-15  # Pa⁻¹

# Facteurs de conversion SUPCRT → SI
const _CAL_TO_J    = 4.184         # 1 cal = 4.184 J
const _BAR_TO_PA   = 1.0e5         # 1 bar = 1e5 Pa
const _ANGST_TO_M  = 1.0e-10       # 1 Å = 1e-10 m
# a₁ [cal/mol/bar] → [J/mol/Pa] : × 4.184 / 1e5
const _A1_CONV     = _CAL_TO_J / _BAR_TO_PA
# a₂ [cal/mol] → [J/mol]
const _A2_CONV     = _CAL_TO_J
# a₃ [cal·K/mol/bar] → [J·K/mol/Pa]
const _A3_CONV     = _CAL_TO_J / _BAR_TO_PA
# a₄ [cal·K/mol] → [J·K/mol]
const _A4_CONV     = _CAL_TO_J
# c₁ [cal/mol/K] → [J/mol/K]
const _C1_CONV     = _CAL_TO_J
# c₂ [cal·K/mol] → [J·K/mol]
const _C2_CONV     = _CAL_TO_J
# ω [cal/mol] → [J/mol]
const _OMG_CONV    = _CAL_TO_J

# ── Structs ───────────────────────────────────────────────────────────────────

"""
    struct HKFParams

Paramètres de l'équation d'état HKF révisée (Tanger & Helgeson 1988) pour une espèce
aqueuse. Stockés en unités SI.

Constructeur : `HKFParams(; a₁, a₂, a₃, a₄, c₁, c₂, ω, units=:supcrt)`

- `units=:supcrt` (défaut) : a₁ [cal/mol/bar], a₂ [cal/mol], a₃ [cal·K/mol/bar],
  a₄ [cal·K/mol], c₁ [cal/mol/K], c₂ [cal·K/mol], ω [cal/mol] — convention SUPCRT92.
- `units=:si` : toutes les valeurs en J et Pa.
"""
struct HKFParams
    a₁::Float64   # J/mol/Pa
    a₂::Float64   # J/mol
    a₃::Float64   # J·K/mol/Pa
    a₄::Float64   # J·K/mol
    c₁::Float64   # J/mol/K
    c₂::Float64   # J·K/mol
    ω::Float64    # J/mol
end

function HKFParams(;
    a₁::Real, a₂::Real, a₃::Real, a₄::Real,
    c₁::Real, c₂::Real, ω::Real,
    units::Symbol = :supcrt,
)
    if units === :supcrt
        return HKFParams(
            Float64(a₁) * _A1_CONV,
            Float64(a₂) * _A2_CONV,
            Float64(a₃) * _A3_CONV,
            Float64(a₄) * _A4_CONV,
            Float64(c₁) * _C1_CONV,
            Float64(c₂) * _C2_CONV,
            Float64(ω)  * _OMG_CONV,
        )
    elseif units === :si
        return HKFParams(
            Float64(a₁), Float64(a₂), Float64(a₃), Float64(a₄),
            Float64(c₁), Float64(c₂), Float64(ω),
        )
    else
        throw(ArgumentError("units doit être :supcrt ou :si, reçu : $units"))
    end
end

"""
    struct HKFReferenceState

État standard de référence à (Tr=298.15 K, Pr=1 bar) pour une espèce aqueuse.
Valeurs en J/mol (ou J/mol/K pour Sr).

Constructeur mot-clé : `HKFReferenceState(Gr=..., Hr=..., Sr=...)`
"""
Base.@kwdef struct HKFReferenceState
    Gr::Float64   # ΔfG⁰(Tr, Pr) [J/mol]
    Hr::Float64   # ΔfH⁰(Tr, Pr) [J/mol]
    Sr::Float64   # S⁰(Tr, Pr)   [J/mol/K]
end

# ── Calcul des propriétés thermodynamiques ────────────────────────────────────

"""
    hkf_Cp(p::HKFParams, T, P; weos=DEFAULT_WATER_EOS) -> Float64

Capacité calorifique standard Cp° [J/mol/K] selon Tanger & Helgeson (1988) :
```
Cp° = c₁ + c₂/(T−Θ)² − ω·T·X
```
"""
function hkf_Cp(p::HKFParams, T::Float64, P::Float64;
                weos::AbstractWaterEOS = DEFAULT_WATER_EOS)::Float64
    bf = born_functions(T, P, weos)
    return p.c₁ + p.c₂ / (T - HKF_Θ)^2 + p.ω * T * bf.X
end

"""
    hkf_V(p::HKFParams, T, P; weos=DEFAULT_WATER_EOS) -> Float64

Volume molaire partiel apparent V° [m³/mol] selon Tanger & Helgeson (1988) :
```
V° = a₁ + a₂/(Ψ+P) + a₃/(T−Θ) + a₄/[(Ψ+P)(T−Θ)] + ω·Q
```
"""
function hkf_V(p::HKFParams, T::Float64, P::Float64;
               weos::AbstractWaterEOS = DEFAULT_WATER_EOS)::Float64
    bf = born_functions(T, P, weos)
    return (
        p.a₁ +
        p.a₂ / (HKF_Ψ + P) +
        p.a₃ / (T - HKF_Θ) +
        p.a₄ / ((HKF_Ψ + P) * (T - HKF_Θ)) +
        p.ω * bf.Q
    )
end

"""
    hkf_S(p::HKFParams, ref::HKFReferenceState, T, P; weos=DEFAULT_WATER_EOS) -> Float64

Entropie standard S° [J/mol/K] selon Tanger & Helgeson (1988).
Intégration de Cp/T de Tr à T (à Pr), puis −∂G/∂T de Pr à P (à T).
"""
function hkf_S(p::HKFParams, ref::HKFReferenceState, T::Float64, P::Float64;
               weos::AbstractWaterEOS = DEFAULT_WATER_EOS)::Float64
    bf_T = born_functions(T, P, weos)
    bf_r = born_functions(HKF_Tr, HKF_Pr, weos)

    Θ = HKF_Θ
    Ψ = HKF_Ψ
    Tr = HKF_Tr
    Pr = HKF_Pr

    # Terme c₁ : intégrale de (c₁/T) dT de Tr à T
    dS_c1 = p.c₁ * log(T / Tr)

    # Terme c₂ : intégrale de c₂/(T·(T−Θ)²) dT de Tr à T
    # = c₂/Θ² × [Θ/(Θ−T) − Θ/(Θ−Tr) + ln(T·(Θ−Tr)/(Tr·(Θ−T)))]
    dS_c2 = (p.c₂ / Θ^2) * (
        Θ / (Θ - T) - Θ / (Θ - Tr) + log(T * (Θ - Tr) / (Tr * (Θ - T)))
    )

    # Terme volumique : dS_vol = −∂G_vol/∂T = +[a₃(P−Pr) + a₄ ln((Ψ+P)/(Ψ+Pr))] / (T−Θ)²
    # (signe positif : G_vol contient [...]/(T−Θ), donc −∂/∂T = +[...]/(T−Θ)²)
    dS_vol =
        (p.a₃ * (P - Pr) + p.a₄ * log((Ψ + P) / (Ψ + Pr))) / (T - Θ)^2

    # Terme Born : ω × ΔY
    dS_born = p.ω * (bf_T.Y - bf_r.Y)

    return ref.Sr + dS_c1 + dS_c2 + dS_vol + dS_born
end

"""
    hkf_H(p::HKFParams, ref::HKFReferenceState, T, P; weos=DEFAULT_WATER_EOS) -> Float64

Enthalpie standard H° [J/mol] selon Tanger & Helgeson (1988).
"""
function hkf_H(p::HKFParams, ref::HKFReferenceState, T::Float64, P::Float64;
               weos::AbstractWaterEOS = DEFAULT_WATER_EOS)::Float64
    bf_T = born_functions(T, P, weos)
    bf_r = born_functions(HKF_Tr, HKF_Pr, weos)
    ε_T = bf_T.ε
    ε_r = bf_r.ε

    Θ = HKF_Θ
    Ψ = HKF_Ψ
    Tr = HKF_Tr
    Pr = HKF_Pr

    # Terme c₁ : ∫Tr→T c₁ dT
    dH_c1 = p.c₁ * (T - Tr)

    # Terme c₂ : ∫Tr→T c₂/(T−Θ)² dT = c₂ × [1/(Θ−T) − 1/(Θ−Tr)]
    dH_c2 = p.c₂ * (1.0 / (Θ - T) - 1.0 / (Θ - Tr))

    # Termes volumiques : ∫Pr→P [V − T(∂V/∂T)] dP'
    # ∫[a₁ + a₂/(Ψ+P') + (a₃+a₄/(Ψ+P')) × (2T−Θ)/(T−Θ)²] dP'
    # (facteur (2T−Θ)/(T−Θ)² provient de 1/(T−Θ) + T/(T−Θ)²)
    Pvol = p.a₃ * (P - Pr) + p.a₄ * log((Ψ + P) / (Ψ + Pr))
    dH_vol = (
        p.a₁ * (P - Pr) +
        p.a₂ * log((Ψ + P) / (Ψ + Pr)) +
        Pvol * (2 * T - Θ) / (T - Θ)^2
    )

    # Terme Born : H_born = ω(1/ε − 1/εr) + ωT·Y(T,P) − ωTr·Yr
    # Dérivé de G_born = ω(1/ε−1/εr) + ωYr(T−Tr) par H = G + T·S
    dH_born = (
        p.ω * (1.0 / ε_T - 1.0 / ε_r) +
        p.ω * T * bf_T.Y -
        p.ω * Tr * bf_r.Y
    )

    return ref.Hr + dH_c1 + dH_c2 + dH_vol + dH_born
end

"""
    hkf_G(p::HKFParams, ref::HKFReferenceState, T, P; weos=DEFAULT_WATER_EOS) -> Float64

Énergie libre de Gibbs standard G° [J/mol] selon Tanger & Helgeson (1988).

Dans les tables SUPCRT, Gr, Hr et Sr sont des quantités indépendantes qui ne satisfont
pas nécessairement Gr = Hr − Tr·Sr (références d'entropie différentes entre espèces et
éléments). La formule correcte est :

    G°(T,P) = H°(T,P) − T·S°(T,P) + (Gr − Hr + Tr·Sr)

Le terme de correction (Gr − Hr + Tr·Sr) est constant et garantit G°(Tr,Pr) = Gr,
tout en préservant les dérivées thermodynamiques ∂G/∂T = −S et ∂G/∂P = V.
"""
function hkf_G(p::HKFParams, ref::HKFReferenceState, T::Float64, P::Float64;
               weos::AbstractWaterEOS = DEFAULT_WATER_EOS)::Float64
    H = hkf_H(p, ref, T, P; weos = weos)
    S = hkf_S(p, ref, T, P; weos = weos)
    # Correction de référence : Gr, Hr, Sr sont indépendants dans les tables SUPCRT
    G_corr = ref.Gr - ref.Hr + HKF_Tr * ref.Sr
    return H - T * S + G_corr
end

# ── Interface callable (compatible ThermoFunction) ────────────────────────────

"""
    struct HKFThermoFunction

Fonction thermodynamique HKF appelable. Interface compatible avec `ThermoFunction` :

```julia
f(; T=298.15, P=1e5, unit=false)
```

`property` peut être `:Cp⁰`, `:V⁰`, `:S⁰`, `:ΔₐH⁰`, `:ΔₐG⁰`.
"""
struct HKFThermoFunction
    params::HKFParams
    refstate::HKFReferenceState
    property::Symbol
    weos::AbstractWaterEOS
end

function (f::HKFThermoFunction)(;
    T::Real = HKF_Tr,
    P::Real = HKF_Pr,
    unit::Bool = false,
)
    T_val = Float64(T)
    P_val = Float64(P)
    p = f.params
    r = f.refstate
    w = f.weos

    val = if f.property === :Cp⁰
        hkf_Cp(p, T_val, P_val; weos = w)
    elseif f.property === :V⁰
        hkf_V(p, T_val, P_val; weos = w)
    elseif f.property === :S⁰
        hkf_S(p, r, T_val, P_val; weos = w)
    elseif f.property === :ΔₐH⁰
        hkf_H(p, r, T_val, P_val; weos = w)
    elseif f.property === :ΔₐG⁰
        hkf_G(p, r, T_val, P_val; weos = w)
    else
        throw(ArgumentError("Propriété inconnue : $(f.property). Valides : :Cp⁰, :V⁰, :S⁰, :ΔₐH⁰, :ΔₐG⁰"))
    end

    return unit ? val * _hkf_unit(f.property) : val
end

"""Unité SI de la propriété thermodynamique HKF."""
function _hkf_unit(prop::Symbol)
    if prop === :Cp⁰ || prop === :S⁰
        return u"J/mol/K"
    elseif prop === :V⁰
        return u"m^3/mol"
    else  # G°, H°
        return u"J/mol"
    end
end

# ── Constructeur de haut niveau ───────────────────────────────────────────────

"""
    build_hkf_functions(params::HKFParams, refstate::HKFReferenceState;
                        weos::AbstractWaterEOS = DEFAULT_WATER_EOS)
    -> OrderedDict{Symbol, HKFThermoFunction}

Construit un dictionnaire de fonctions thermodynamiques HKF pour une espèce aqueuse.
Interface identique à `build_thermo_functions` dans `thermo_models.jl`.

Clés : `:Cp⁰`, `:V⁰`, `:S⁰`, `:ΔₐH⁰`, `:ΔₐG⁰`

```julia
params = HKFParams(; a₁=..., a₂=..., a₃=..., a₄=..., c₁=..., c₂=..., ω=..., units=:supcrt)
ref    = HKFReferenceState(Gr=-261.9e3, Hr=-240.3e3, Sr=59.0)
fns    = build_hkf_functions(params, ref)
fns[:ΔₐG⁰](T=298.15, P=1e5)   # → -261.9e3 J/mol
```
"""
function build_hkf_functions(
    params::HKFParams,
    refstate::HKFReferenceState;
    weos::AbstractWaterEOS = DEFAULT_WATER_EOS,
)
    props = (:Cp⁰, :V⁰, :S⁰, :ΔₐH⁰, :ΔₐG⁰)
    return OrderedDict{Symbol, HKFThermoFunction}(
        prop => HKFThermoFunction(params, refstate, prop, weos) for prop in props
    )
end
