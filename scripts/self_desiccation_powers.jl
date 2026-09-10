# =============================================================================
# self_desiccation_powers.jl
#
# Where does Powers' 0.42 come from?
#
# Powers' rule of thumb caps the degree of hydration of a sealed paste at
# `α_max = (w/c) / 0.42`. This script takes that coefficient apart into a water
# budget, supplies each term from an independent source — the chemistry from a
# Gibbs energy minimization, the arrest saturation from a published desorption
# isotherm — and reports what the two together imply.
#
# Nothing here is fitted to 0.42.
#
# Usage:
#   julia --project=scripts scripts/self_desiccation_powers.jl
#   or from the REPL:  include("scripts/self_desiccation_powers.jl")
#
# Companion to the tutorial "Self-desiccation and Powers' 0.42".
# =============================================================================

import Pkg
Pkg.activate(@__DIR__; io = devnull)

using ChemistryLab
using DynamicQuantities
using OptimaSolver
using Printf

# ── The cement, and the paste ────────────────────────────────────────────────
#
# Baroghel-Bouny et al. (1999), Table 2: the cement whose desorption isotherm
# this script uses further down. Taking the isotherm from one paper and the
# clinker from another would compare two materials.

const COMPO = [
    "C3S" => 0.5728, "C2S" => 0.2398, "C3A" => 0.0303,
    "C4AF" => 0.0759, "Gp" => 0.0439, "Cal" => 0.0184,
]
const CMASS = sum(last.(COMPO))      # 0.9811 — the rest is free lime and alkalis
const WC = 0.34                      # their mix CO
const M_H2O = 0.0180153              # kg/mol
const RHO_W = 1.0e3                  # kg/m³, for the water volume

const PHASES = split(
    "C3S C2S C3A C4AF Gp Anh Cal Portlandite Jennite H2O@ ettringite " *
        "monosulphate12 C3AH6 C3FH6 C4FH13 monocarbonate"
)

substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
species = speciation(substances, PHASES; aggregate_state = [AS_AQUEOUS])
cs = ChemicalSystem(species, CEMDATA_PRIMARIES)

const IW = only(cs.idx_solvent)
const IC = [findfirst(s -> symbol(s) == sym, cs.species) for (sym, _) in COMPO]

"""
    paste(α) -> ChemicalState

A paste in which a fraction `α` of the cement has been made available to react,
with **all** of the mixing water. The unreacted `1 - α` is left out of the
minimization and added back afterwards for the volume balance: an equilibrium
calculation has no notion of a reaction that has not happened, so the degree of
hydration has to be imposed. This is the construction of Lothenbach & Winnefeld
(2006).
"""
function paste(α)
    mtot = CMASS + WC * CMASS
    st = ChemicalState(cs)
    for (sym, mfrac) in COMPO
        set_quantity!(st, sym, α * mfrac / mtot * u"kg")
    end
    set_quantity!(st, "H2O@", WC * CMASS / mtot * u"kg")
    V = volume(st)
    set_quantity!(st, "H+", 1.0e-7u"mol/L" * V.liquid)
    set_quantity!(st, "OH-", 1.0e-7u"mol/L" * V.liquid)
    return st
end

cement_mass(st) = sum(ustrip(us"kg", st.n[i] * cs.species[i][:M]) for i in IC)

"""
    budget(α) -> NamedTuple

The water budget at imposed degree of hydration `α`, per gram of cement:

  - `b`   — water bound in the hydrate formulae, g/g, divided by `α`
  - `s`   — chemical shrinkage (the empty porosity), cm³/g, divided by `α`
  - `w_free` — water left as pore solution, g/g
  - `porosity`, `certified`

Dividing by `α` is the test of the linear ansatz the budget below assumes: if `b`
and `s` come out independent of `α`, the assumption holds for this system.
"""
function budget(α)
    fresh = paste(1.0)
    mc = cement_mass(fresh)
    w_tot = ustrip(us"mol", fresh.n[IW]) * M_H2O
    eq, cert = equilibrate_certified(paste(α))

    n = collect(eq.n)
    for i in IC
        n[i] += (1 - α) * fresh.n[i]
    end
    final = ChemicalState(cs, n)
    ϕ = porosity(final, fresh)

    w_free = ustrip(us"mol", eq.n[IW]) * M_H2O
    V_ref = ustrip(us"m^3", volume(fresh).total)
    return (;
        b = (w_tot - w_free) / mc / α,
        s = ϕ.void * V_ref / mc / α * 1.0e3,     # m³/kg → cm³/g
        w_free = w_free / mc,
        porosity = ϕ.total,
        certified = cert.optimal,
    )
end

# ── 1. The chemistry: b and s ────────────────────────────────────────────────

println("="^78)
println("1. THE WATER BUDGET, FROM THE CHEMISTRY")
println("="^78)
println("  alpha   b (g/g)   s (cm3/g)   free water (g/g)   porosity   certified")
for α in (0.55, 0.6, 0.65, 0.7, 0.8)
    r = budget(α)
    @printf(
        "%7.2f  %9.4f  %10.4f  %16.4f  %9.4f   %s\n",
        α, r.b, r.s, r.w_free, r.porosity, r.certified
    )
end

const REF = budget(0.65)
const B_MODEL = REF.b
const S_SHRINK = REF.s
@printf(
    "\nb = %.4f g/g (formula water)   s = %.4f cm3/g (chemical shrinkage)\n",
    B_MODEL, S_SHRINK
)

# ── 2. The retention curve, from measurement ─────────────────────────────────
#
# Baroghel-Bouny et al. (1999), their Eq. (20) and Table 5, mix CO. They write
# the van Genuchten expression with `b = 1/m`, so their `b = 2.1684` enters here
# as `m = 1/2.1684`. Getting that inversion wrong is silent.

const CO_CURVE = VanGenuchten(; a = 37.5479e6, m = 1 / 2.1684)
const V_M_WATER = 1.807e-5      # m³/mol
const T_K = 298.15

water_activity_at(S) = water_activity(CO_CURVE, S; V_m = V_M_WATER, T = T_K)

"""
    saturation_at(rh) -> Float64

The degree of saturation at which the measured isotherm holds its water at
activity `rh`. Bisection, because the van Genuchten form is not invertible in
closed form for this exponent pair.
"""
function saturation_at(rh)
    lo, hi = 1.0e-4, 1.0 - 1.0e-12
    for _ in 1:60
        mid = (lo + hi) / 2
        water_activity_at(mid) > rh ? (hi = mid) : (lo = mid)
    end
    return (lo + hi) / 2
end

println("\n", "="^78)
println("2. THE MEASURED DESORPTION ISOTHERM (Baroghel-Bouny et al. 1999, mix CO)")
println("="^78)
println("    RH     S*      capillary pressure (MPa)   Kelvin radius (nm)")
for rh in (0.75, 0.8, 0.85, 0.9, 0.95)
    S = saturation_at(rh)
    @printf(
        "  %5.2f  %6.4f  %24.2f  %18.2f\n", rh, S,
        capillary_pressure(CO_CURVE, S) / 1.0e6,
        kelvin_radius(rh; γ = 0.0728u"N/m", V_m = 1.807e-5u"m^3/mol", T = 298.15u"K") * 1.0e9
    )
end

# ── 3. The budget closed: Powers' coefficient ────────────────────────────────
#
#   per gram of cement, at degree of hydration α:
#       bound water        = b α                       (g)
#       free water         = w/c − b α                 (g)
#       liquid volume      = (w/c − b α) / ρ_w         (cm³)
#       empty volume       = s α                       (cm³)   ← Le Chatelier
#       S                  = V_liq / (V_liq + V_void)
#
#   arrest at S = S* gives  α_max = (w/c) / k  with  k = b + s S*/(1 − S*)

powers_k(b, s, Sstar) = b + s * Sstar / (1 - Sstar)

function invert_k(b, s, k_target)
    lo, hi = 1.0e-6, 1.0 - 1.0e-12
    for _ in 1:80
        mid = (lo + hi) / 2
        powers_k(b, s, mid) > k_target ? (hi = mid) : (lo = mid)
    end
    return (lo + hi) / 2
end

const S_AT_80 = saturation_at(0.8)
const B_POWERS = 0.23     # Powers' non-evaporable water, g/g — his own number

println("\n", "="^78)
println("3. THE BUDGET CLOSED")
println("="^78)
@printf("S* at RH 0.80 = %.4f,   so S*/(1-S*) = %.4f\n\n", S_AT_80, S_AT_80 / (1 - S_AT_80))

@printf("FORWARD, with the model's formula water b = %.4f:\n", B_MODEL)
@printf("    k = %.4f          (Powers: 0.42)\n", powers_k(B_MODEL, S_SHRINK, S_AT_80))
@printf(
    "    alpha_max(w/c = %.2f) = %.4f   (Powers: %.4f)\n\n",
    WC, WC / powers_k(B_MODEL, S_SHRINK, S_AT_80), min(1.0, WC / 0.42)
)

@printf("FORWARD, with Powers' own non-evaporable water b = %.2f:\n", B_POWERS)
@printf("    k = %.4f          (Powers: 0.42)\n", powers_k(B_POWERS, S_SHRINK, S_AT_80))
@printf(
    "    alpha_max(w/c = %.2f) = %.4f   (Powers: %.4f)\n\n",
    WC, min(1.0, WC / powers_k(B_POWERS, S_SHRINK, S_AT_80)), min(1.0, WC / 0.42)
)

println("INVERTED — what internal humidity does Powers' 0.42 imply?")
for (label, b) in (("the model's formula water", B_MODEL), ("Powers' own w_n", B_POWERS))
    S = invert_k(b, S_SHRINK, 0.42)
    @printf(
        "    b = %.4f (%-26s):  S* = %.4f  =>  RH = %.4f\n",
        b, label, S, water_activity_at(S)
    )
end

# ── 4. The discrepancy, attributed ───────────────────────────────────────────

println("\n", "="^78)
println("4. THE DISCREPANCY")
println("="^78)
@printf("model formula water  %.4f g/g\n", B_MODEL)
@printf("Powers' w_n          %.4f g/g   (defined by D-drying)\n", B_POWERS)
@printf(
    "difference           %.4f g/g   — interlayer water CEMDATA18 writes\n",
    B_MODEL - B_POWERS
)
println("                                  into the C-S-H formula and D-drying removes")
@printf("\nPowers' own split:   w_n %.2f + gel %.2f = %.2f\n", 0.23, 0.19, 0.42)
@printf(
    "this budget:         b   %.4f + s S*/(1-S*) %.4f = %.4f\n",
    B_MODEL, S_SHRINK * S_AT_80 / (1 - S_AT_80), powers_k(B_MODEL, S_SHRINK, S_AT_80)
)

@printf(
    "\nSensitivity: dk/dS* = s/(1-S*)^2 = %.3f at S* = %.4f\n",
    S_SHRINK / (1 - S_AT_80)^2, S_AT_80
)
println("  so a tenth of a point of saturation is worth 0.14 in k, and the")
println("  arrest humidity has to be known to about a point to pin k to 0.01.")

println("\nMeasured porosity of their mix CO (their Table 4): 30.3 %")
@printf("This model at alpha = 0.65: %.1f %%\n", 100 * REF.porosity)
