# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities
using LinearAlgebra

# ── Abstract activity model ───────────────────────────────────────────────────

"""
    abstract type AbstractActivityModel end

Base type for all activity models: the correction that turns a concentration
into a chemical potential.

## What an activity model is for

The chemical potential of a species is

```math
\\mu_i = \\mu_i^\\circ + RT \\ln a_i ,
```

and the whole of the equilibrium calculation — [`build_potentials`](@ref), the
Gibbs minimization, the certificate — depends on nothing else about the
solution. So an activity model is the single place where a real solution stops
being ideal, and it has to answer **two** questions, not one:

1. the **solute** activity, `a_i = γ_i m_i` (or `γ_i c_i`), where `γ_i` corrects
   for the fact that an ion in a solution of other ions is not in the same state
   as one alone at the same concentration;
2. the **solvent** activity `a_w`, which is not a small correction on the way to
   the same answer. It is what a hydrate reaction consumes, and in a cement
   paste with little mixing water it is the quantity that decides how far the
   reaction goes.

The second is the one to look at when choosing a model here, because the three
built-in models treat it in three different ways and only one of them is
defensible past a dilute solution. See
[Activity models](@ref sec-activity-models) for the comparison.

## The interface

A concrete subtype must implement

```julia
activity_model(cs::ChemicalSystem, model::YourModel) -> lna    # lna(n, p)
concentration_scale(model::YourModel)                          # :molality | :molarity
```

`lna(n, p)` receives the **full mole vector** `n`, indexed like `cs.species`,
and the parameter tuple `p` carrying at least `ϵ` and usually `T`, `P` and
`ΔₐG⁰overT`; it returns `ln aᵢ` for every species — solutes, solvent, pure
crystals (`0`), gases, and solid-solution end-members. Three properties are
required rather than nice to have:

  - it is differentiated by `ForwardDiff` at every Newton step, so the output
    element type must follow `n` and any regularization must be smooth;
  - it must call `_solid_solution_lna!`, or the end-members of a solid solution
    silently get `ln a = 0`;
  - `concentration_scale` has no fallback: without it the aqueous accessors
    raise a `MethodError` rather than guessing a convention.

One discipline is worth knowing before adding a model: the scalar formula lives
in `_log10γ_ion` and `_log10γ_neutral`, which the closure calls *and*
[`activity_coefficients`](@ref) calls, so a model states its algebra once and the
accessor cannot drift from what the solver used. A model whose coefficients are
not a per-species function of `(z, å, I)` — an ion-interaction model, where the
sum runs over pairs — does not fit that pair of helpers and supplies its own
vector-valued path instead.

See also: [`DiluteSolutionModel`](@ref), [`HKFActivityModel`](@ref),
[`DaviesActivityModel`](@ref), [`concentration_scale`](@ref).
"""
abstract type AbstractActivityModel end

# ── Concrete models ───────────────────────────────────────────────────────────

"""
    struct DiluteSolutionModel <: AbstractActivityModel

The ideal dilute solution: every activity coefficient is exactly 1, and an
activity is a concentration.

```math
\\gamma_i \\equiv 1 , \\qquad
a_i = \\frac{c_i}{c^\\circ}\\ (c^\\circ = 1\\ \\mathrm{mol/L}) , \\qquad
a_w = x_w
```

| phase | law | expression |
|:--|:--|:--|
| solvent | Raoult | `ln a = ln x_w`, the **mole fraction** |
| aqueous solutes | Henry | `ln a = ln(cᵢ/c°)`, `c° = 1 mol/L` |
| pure crystals | — | `ln a = 0` |
| gas | ideal mixture | `ln a = ln xᵢ` |
| solid-solution end-members | ideal mixing | `ln a = ln xᵢ` within the phase |

## The physics, and where it runs out

Ideality means an ion does not notice the others. That is true in the limit of
infinite dilution and stops being true as soon as the ionic cloud around an ion
has a measurable effect on its energy — that is, from a few millimolal upwards
for a charged species. Take `I ≲ 0.01 mol/kg` as the range in which the answer
is the answer, and treat anything above as a screening calculation.

!!! warning "Its water activity is a mole fraction, and that is the real limit"
    `a_w = x_w` is Raoult's law, which counts molecules and knows nothing about
    what they are. A cement pore solution at `x_w = 0.99` gets `a_w = 0.99`
    whatever it holds in solution, and a paste short of mixing water gets an
    `a_w` that follows the *amount* of water and not its state. If the water
    activity matters for what is being computed — and in a hydrating paste it
    decides how far the reaction goes — this model is the wrong one, and
    [`HKFActivityModel`](@ref) is the least that will do.

It is nonetheless the **default**, for two reasons that are about the solve and
not about the chemistry: it is exact in the dilute limit, and it makes the
log-activity linear in `ln n`, which is the best-conditioned objective the
minimizer will ever see. Start here, then change the model and see whether the
answer moves.

See also: [`HKFActivityModel`](@ref), [`DaviesActivityModel`](@ref),
[Activity models](@ref sec-activity-models).
"""
struct DiluteSolutionModel <: AbstractActivityModel end

# ── Activity model factory ────────────────────────────────────────────────────

"""
    activity_model(cs::ChemicalSystem, ::DiluteSolutionModel) -> Function

Return a closure `lna(n, p) -> Vector{Float64}` computing the vector of
log-activities for the dilute ideal solution model.

The returned function has signature `lna(n, p)` where:
- `n`: dimensionless mole vector (same indexing as `cs.species`)
- `p`: `NamedTuple` containing at least `ϵ` (floor value to avoid log(0))

Solid-solution end-members (class `SC_SSENDMEMBER`) receive `ln aᵢ = ln xᵢ`
where `xᵢ = nᵢ / Σnⱼ` within the same solid-solution phase.

All quantities are dimensionless — units are stripped at construction time.
"""
function activity_model(cs::ChemicalSystem, ::DiluteSolutionModel)

    has_aqueous = !isempty(cs.idx_solvent)
    idx_solvent = has_aqueous ? only(cs.idx_solvent) : 0
    idx_solutes = cs.idx_solutes
    idx_gas = cs.idx_gas

    ln_c_solvent = if has_aqueous
        M_solvent = ustrip(us"kg/mol", cs.species[idx_solvent][:M])
        log(1.0 / M_solvent)    # c° = ρ/M ≈ 1/M (ρ ≈ 1 kg/L)
    else
        0.0
    end

    ss_groups = cs.ss_groups
    has_ss = !isempty(ss_groups)
    has_gas = !isempty(idx_gas)
    ss_models = has_ss ? map(ss -> ss.model, cs.solid_solutions) : nothing

    function lna(n::AbstractVector, p)
        ϵ = p.ϵ
        _n = max.(n, ϵ)     # ϵ::Float64 — promotion vers Dual automatique si n est Dual

        out = zeros(eltype(_n), length(_n))

        if has_aqueous
            # n_aqueous ≥ ϵ > 0 always (because _n[i] ≥ ϵ), so no iszero guard needed
            n_aqueous = _n[idx_solvent] + sum((_n[i] for i in idx_solutes); init = zero(eltype(_n)))
            out[idx_solvent] = log(_n[idx_solvent] / n_aqueous)
            @inbounds for i in idx_solutes
                out[i] = log(_n[i] / _n[idx_solvent]) + ln_c_solvent
            end
        end

        if has_gas
            n_gas = sum((_n[i] for i in idx_gas); init = zero(eltype(_n)))
            @inbounds for i in idx_gas
                out[i] = log(_n[i] / n_gas)
            end
        end

        if has_ss
            T_val = hasproperty(p, :T) ? p.T : 298.15
            _solid_solution_lna!(out, _n, ss_groups, ss_models, T_val, ϵ)
        end

        return out
    end

    return lna
end

# ── Potential builder ─────────────────────────────────────────────────────────

"""
    build_potentials(cs::ChemicalSystem, model::AbstractActivityModel) -> Function

Return a closure `μ(n, p) -> Vector{Float64}` computing dimensionless chemical
potentials `μ_i / RT` for all species.

``\\mu_i / RT = \\Delta_a G_i^0 / RT + \\ln a_i``

The returned function is compatible with SciML solvers:
- `n`: dimensionless mole vector
- `p`: `NamedTuple` containing:
  - `ΔₐG⁰overT`: vector of standard Gibbs energies of formation divided by RT
  - `ϵ`: regularization floor (e.g. `1e-30`)

All quantities are dimensionless — caller is responsible for stripping units
from `ΔₐG⁰overT` before passing them in `p`.

# Examples
```jldoctest
julia> cs = ChemicalSystem([
           Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("Na+"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
       ]);

julia> μ = build_potentials(cs, DiluteSolutionModel());

julia> n = [55.5, 0.1];

julia> p = (ΔₐG⁰overT = [-95.6, -105.6], ϵ = 1e-30);

julia> length(μ(n, p)) == 2
true
```
"""
function build_potentials(cs::ChemicalSystem, model::AbstractActivityModel)

    # Build the activity closure once — captures precomputed indices and constants
    lna = activity_model(cs, model)

    function μ(n::AbstractVector, p)
        return p.ΔₐG⁰overT .+ lna(n, p)       # μ_i/RT = ΔₐG⁰_i/RT + ln(a_i)
    end

    return μ
end

# ── HKF / Debye-Hückel B-dot activity model ──────────────────────────────────

"""
    REJ_HKF::Dict{String,Float64}

Effective electrostatic radii åᵢ [Å] for aqueous ions from
Helgeson, Kirkham & Flowers (1981), *Am. J. Sci.* **281**, Table 3.

Keys are PHREEQC-format formula strings (e.g. `"Na+"`, `"Ca+2"`, `"SO4-2"`).
Used by [`HKFActivityModel`](@ref) with priority 2 in the radius lookup chain:
`sp[:å]` > `REJ_HKF` > [`REJ_CHARGE_DEFAULT`](@ref) > `model.å_default`.

See also: [`REJ_CHARGE_DEFAULT`](@ref), [`HKFActivityModel`](@ref).
"""
const REJ_HKF = Dict{String, Float64}(
    "H+" => 3.08, "Li+" => 1.64, "Na+" => 1.91, "K+" => 2.27,
    "Rb+" => 2.41, "Cs+" => 2.61, "NH4+" => 2.31, "Ag+" => 2.2,
    "Mg+2" => 2.54, "Ca+2" => 2.87, "Sr+2" => 3.0, "Ba+2" => 3.22,
    "Fe+2" => 2.62, "Al+3" => 3.33, "Fe+3" => 3.46, "La+3" => 3.96,
    "F-" => 1.33, "Cl-" => 1.81, "Br-" => 1.96, "I-" => 2.2,
    "OH-" => 1.4, "HS-" => 1.84, "NO3-" => 2.81, "HCO3-" => 2.1,
    "HSO4-" => 2.37, "SO4-2" => 3.15, "CO3-2" => 2.81,
)

"""
    REJ_CHARGE_DEFAULT::Dict{Int,Float64}

Fallback effective electrostatic radii åᵢ [Å] indexed by formal charge,
from ToughReact V2 (Xu et al. 2011, Table A2; after Helgeson et al. 1981).

Used by [`HKFActivityModel`](@ref) with priority 3 in the radius lookup chain:
`sp[:å]` > [`REJ_HKF`](@ref) > `REJ_CHARGE_DEFAULT` > `model.å_default`.

See also: [`REJ_HKF`](@ref), [`HKFActivityModel`](@ref).
"""
const REJ_CHARGE_DEFAULT = Dict{Int, Float64}(
    -3 => 4.2, -2 => 3.0, -1 => 1.81,
    1 => 2.31, 2 => 2.8, 3 => 3.6, 4 => 4.5,
)

# ── Internal helpers ──────────────────────────────────────────────────────────

"""
    _hkf_sigma(x) -> Real

Compute the σ function used in the osmotic coefficient formula
(Helgeson et al. 1981, Eq. 132–137):

```
σ(x) = (3/x³)(x − 2 ln(1+x) − 1/(1+x) + 1)
```

For `|x| < 1e-3`, a Taylor series `1 − (3/2)x + (9/5)x²` is used to avoid
catastrophic cancellation. Both branches agree to O(x³) at the threshold,
so the gradient is continuous.

AD-compatible: branching is on `ForwardDiff.value(x)`, not on the Dual itself.
"""
function _hkf_sigma(x::T) where {T <: Real}
    if abs(_primal(x)) < 1.0e-3
        return one(T) - (3 // 2) * x + (9 // 5) * x^2
    else
        return (3 / x^3) * (x - 2 * log1p(x) - one(T) / (one(T) + x) + one(T))
    end
end

"""
    hkf_debye_huckel_params(T_K, P_Pa) -> NamedTuple{(:A, :B)}

Compute the Debye-Hückel A and B parameters from the water density ρ [g/cm³]
and dielectric constant εᵣ at temperature `T_K` (K) and pressure `P_Pa` (Pa).

Formulas (Helgeson et al. 1981):
```
A(T,P) = 1.824829238×10⁶ × √ρ / (εᵣ T)^(3/2)    [(kg/mol)^(1/2)]
B(T,P) = 50.29158649      × √ρ / √(εᵣ T)          [Å⁻¹ (kg/mol)^(1/2)]
```
where ρ is in g/cm³.

Uses `water_thermo_props` (HGK equation of state) and
`water_electro_props_jn` (Johnson-Norton dielectric constant).

AD-compatible (ForwardDiff-safe). Returns a `NamedTuple` `(A=..., B=...)`.

# Examples
```jldoctest
julia> p = hkf_debye_huckel_params(298.15, 1e5);

julia> isapprox(p.A, 0.5114; rtol=1e-3)
true

julia> isapprox(p.B, 0.3288; rtol=1e-3)
true
```
"""
function hkf_debye_huckel_params(T_K, P_Pa)
    T_K, P_Pa = promote(T_K, P_Pa)
    wtp = water_thermo_props(T_K, P_Pa)
    wep = water_electro_props_jn(T_K, P_Pa, wtp)
    ρ_gcm3 = wtp.D / 1000                       # kg/m³ → g/cm³
    εT = wep.epsilon * T_K                   # dimensionless × K
    A = 1.824829238e6 * sqrt(ρ_gcm3) / εT^(3 // 2)
    B = 50.29158649 * sqrt(ρ_gcm3) / sqrt(εT)
    return (A = A, B = B)
end

# Internal: Setschenow (salting-out) coefficient lookup for a neutral aqueous
# species. `model.Kₙ` is one coefficient for every neutral species, which is what
# the B-dot literature assumes; a per-species value is what actually varies —
# CO2(aq), the noble gases and the neutral silicates are not equally salted out.
# Priority: `sp[:Kₙ]`, then the model's global value. No table is shipped,
# because a table of Setschenow coefficients is data, and data belongs in a
# database or in the caller's hands, not hard-coded here.
function _setschenow(sp::AbstractSpecies, model)
    haskey(properties(sp), :Kₙ) && return float(sp[:Kₙ])
    return float(model.Kₙ)
end

# Internal: ionic radius priority lookup. A model-level `å` short-circuits the
# whole chain — that is the point of it: a common radius must not be silently
# overridden by a per-species table entry.
function _hkf_lookup_å(sp::AbstractSpecies, model)
    if hasproperty(model, :å) && model.å !== nothing
        return float(model.å)
    end
    if haskey(properties(sp), :å)
        v = sp[:å]
        return v isa Number ? float(v) : float(safe_ustrip(1.0u"Å", v))
    end
    pf = phreeqc(formula(sp))
    haskey(REJ_HKF, pf)           && return REJ_HKF[pf]
    z = Int(charge(sp))
    haskey(REJ_CHARGE_DEFAULT, z) && return REJ_CHARGE_DEFAULT[z]
    return float(model.å_default)
end

# ── HKFActivityModel ──────────────────────────────────────────────────────────

"""
    struct HKFActivityModel{T<:Real} <: AbstractActivityModel

The extended Debye-Hückel model with a B-dot term — the workhorse of aqueous
geochemistry, and the model PHREEQC and EQ3/6 use
([Helgeson1969](@cite), [Helgeson1981](@cite), [ParkhurstAppelo2013](@cite)).

# Formulas

Ionic species (`z ≠ 0`), on the molality scale with `m° = 1 mol/kg`:

```math
\\log_{10}\\gamma_i =
   -\\,\\frac{A\\,z_i^2\\sqrt{I}}{1 + B\\,\\mathring{a}_i\\sqrt{I}}
   \\;+\\; \\dot{B}\\,I ,
\\qquad
I = \\tfrac{1}{2}\\sum_j m_j z_j^2
```

Neutral aqueous species (`z = 0`), the Setschenow form:

```math
\\log_{10}\\gamma_i = K_n I
```

and in both cases `\\ln a_i = \\ln 10 \\cdot \\log_{10}\\gamma_i + \\ln m_i`.

The **water** activity comes from the osmotic coefficient `φ`, so that it is
consistent with the same parameters rather than being a separate assumption:

```math
\\ln a_w = -M_w \\varphi \\sum_j m_j ,
\\qquad
\\varphi = 1 - \\frac{A\\ln 10}{3}\\,\\frac{\\sum_j m_j z_j^2}{\\sum_j m_j}\\,
              \\sqrt{I}\\,\\sigma\\!\\left(B\\,\\mathring{a}_{\\mathrm{eff}}\\sqrt{I}\\right)
            + \\frac{\\dot{B}\\ln 10}{2} I
```

with `σ(x) = (3/x³)(x − 2\\ln(1+x) − 1/(1+x) + 1)` ([`_hkf_sigma`](@ref),
Helgeson et al. 1981 Eqs. 132–137).

# The physics each term carries

  - **`−A z² √I`** — an ion polarizes the solution around it, and the resulting
    ionic cloud screens its charge and lowers its energy. Hence `γ < 1`: an ion
    in an electrolyte is *more* stable than one alone at the same molality. The
    `√I` dependence is not fitted, it is the Debye-Hückel limiting law, exact as
    `I → 0`, and `A` is fixed by the dielectric constant and density of water.
  - **`1 + B å √I`** — the cloud cannot approach closer than the ion's own
    size, which cuts the screening off. This is the one term that carries
    something about the *identity* of the ion (through `å`), and it is what
    extends a law valid at millimolal to roughly `1 mol/kg`.
  - **`+ Ḃ I`** — empirical, linear, positive, and the reason `γ` turns back
    upwards at high ionic strength. It is not a physical term that was derived
    and then measured: [AndersonCrerar1993](@cite) (§17.7.1, pp. 445-446)
    describe it as a **deviation function**, defined by Helgeson as the
    difference between the observed activity coefficient of an electrolyte —
    NaCl — and what the Debye-Hückel expression predicts for it. So it carries
    short-range ion-solvent and ion-ion interaction *and* whatever the first two
    terms failed to capture, together, in one number fitted to one salt. That is
    why this model has a ceiling rather than an asymptote, and why the ceiling
    is somewhere around a molal rather than at a sharp value.
  - **`K_n I`** for a neutral species — salting out. Water engaged around ions
    is water unavailable to solvate a neutral molecule, so its activity rises
    with `I` and its solubility falls. `CO₂(aq)` is the case that matters here.

# Inputs, their defaults, and where each default comes from

| field | default | unit | provenance |
|:--|:--|:--|:--|
| `A` | 0.5114 | (kg/mol)^½ | [Helgeson1981](@cite) Table 1 at 25 °C / 1 bar — **and** reproduced to `1e-3` by [`hkf_debye_huckel_params`](@ref) from this package's own water model, which is the check in `test/activities.jl` |
| `B` | 0.3288 | Å⁻¹(kg/mol)^½ | same |
| `Ḃ` | 0.041 | kg/mol | the value conventionally carried for a NaCl-dominated solution at 25 °C. What the term *is* has a source ([AndersonCrerar1993](@cite) §17.7.1); **this particular number has none recorded in this package**, so treat it as a convention rather than a measurement |
| `Kₙ` | 0.1 | kg/mol | a generic salting-out coefficient, **no source recorded**; overridden per species by `sp[:Kₙ]`, which is how `CO₂(aq)` gets its own |
| `å_default` | 3.72 | Å | last resort, reached only for a charge no table covers (`|z| ≥ 5`). **No source recorded** |
| `å` | `nothing` | Å | one common radius for every ion, overriding the tables. `å = 0` collapses the denominator and gives the limiting law plus `Ḃ I` |
| `temperature_dependent` | `false` | — | recompute `A` and `B` from `p.T`, `p.P` at every call (needs `T` and `P` in `p`) |

Per-ion radii come from [`REJ_HKF`](@ref) ([Helgeson1981](@cite) Table 3) and,
failing that, from [`REJ_CHARGE_DEFAULT`](@ref) ([Xu2011](@cite) Table A2).

!!! note "Three of these defaults are conventions, not data"
    `Ḃ`, `Kₙ` and `å_default` are numbers this package carries without a source
    to point at. They are in the range everyone uses and they are almost
    certainly right, but the honest statement is that they have not been traced,
    and each is a keyword away from being replaced. `A`, `B` and the radius
    tables are traceable, and `A` and `B` are additionally *derived* here rather
    than tabulated.

# What is approximated in the water activity

The `γᵢ` use a per-ion radius; the osmotic coefficient uses **one**
charge-weighted mean radius `å_eff = Σ mᵢzᵢ²åᵢ / Σ mᵢzᵢ²`. So `a_w` is not
exactly the Gibbs-Duhem integral of the `γᵢ` this model returns when the ions
differ in size, and the inconsistency is measurable: the Gibbs-Duhem residual on
0.3 mol/kg NaCl is a few parts in a thousand (`test/activities.jl` asserts
`< 5e-3`), where an exactly consistent model would sit at solver tolerance.
Å-level differences between Na⁺ (1.91) and Cl⁻ (1.81) are enough to produce it.

That is acceptable for a pore solution at 0.1–0.5 mol/kg and it is the first
thing that breaks as the solution concentrates. It is also *structural*: no
choice of `Ḃ` repairs it, because the defect is in using a mean radius at all.

# Ionic radius lookup

The effective radius `åᵢ` is resolved in order:

1. `model.å` — a common radius for every ion, when given. Short-circuits the
   rest of the chain, so a per-species table entry cannot silently override it.
2. `sp[:å]` — explicit value in the species properties dict.
3. [`REJ_HKF`](@ref) — Helgeson et al. (1981) Table 3, keyed by PHREEQC formula.
4. [`REJ_CHARGE_DEFAULT`](@ref) — fallback by formal charge.
5. `model.å_default` — reached only for a charge no table covers, i.e. `|z| ≥ 5`.
   It is **not** a way to impose a common radius; pass `å` for that.

# Valid range

`I ≲ 1 mol/kg` for the activity coefficients, and less than that for the water
activity, for the reason above. Beyond a few mol/kg the `Ḃ I` term is doing work
it was never fitted to do, and nothing in the formula announces it: the model
goes on returning finite, plausible numbers. A cement pore solution normally
sits at 0.1–0.5 mol/kg and is comfortably inside; a paste short of mixing water
is not — see [`solvent_fraction`](@ref) and
[Activity models](@ref sec-activity-models).

# Examples
```julia
# Default model at 25 °C / 1 bar (fixed A, B)
model = HKFActivityModel()

# Temperature-dependent A and B (recomputed at each solve)
model_tdep = HKFActivityModel(temperature_dependent=true)

state_eq = equilibrate(state; model=HKFActivityModel())
```

See also: [`DaviesActivityModel`](@ref), [`hkf_debye_huckel_params`](@ref),
[`activity_coefficients`](@ref), [`ionic_strength`](@ref).
"""
struct HKFActivityModel{T <: Real} <: AbstractActivityModel
    A::T
    B::T
    Ḃ::T
    Kₙ::T
    å_default::T
    å::Union{Nothing, T}
    temperature_dependent::Bool
end

"""
    HKFActivityModel(; A=0.5114, B=0.3288, Ḃ=0.041, Kₙ=0.1, å_default=3.72,
                       å=nothing, temperature_dependent=false) -> HKFActivityModel

Construct an [`HKFActivityModel`](@ref) with the given parameters.

Default values are from Helgeson et al. (1981), Table 1, at 25 °C / 1 bar.

`å` imposes **one common** effective radius on every charged aqueous species,
overriding the per-species tables. Use it to reproduce a published model that
was run with a single ion-size parameter — which is what GEM-Selektor, PHREEQC's
`-gamma` and most cement models do. Note that `å_default` does **not** do this:
it is only the last resort of the lookup chain and is reached only for charges
no table covers. `å = 0` gives the Debye-Hückel limiting law plus the B-dot
term.

# Examples

```julia
# The package default: per-species radii from REJ_HKF, EQ3/6 NaCl B-dot.
HKFActivityModel()

# One common radius of 3.72 Å for every ion.
HKFActivityModel(å = 3.72)

# The Debye-Hückel limiting law with a KOH-background B-dot, which is what a
# GEM-Selektor CEMDATA18 run of a Portland cement uses: CEMDATA18 carries no
# ion-size parameter, so GEMS starts from å = 0, and the B-dot term is not
# applied to neutral species.
HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0)
```
"""
function HKFActivityModel(;
        A::Real = 0.5114,
        B::Real = 0.3288,
        Ḃ::Real = 0.041,
        Kₙ::Real = 0.1,
        å_default::Real = 3.72,
        å::Union{Nothing, Real} = nothing,
        temperature_dependent::Bool = false,
    )
    if å === nothing
        vals = promote(A, B, Ḃ, Kₙ, å_default)
        T = eltype(vals)
        return HKFActivityModel{T}(vals..., nothing, temperature_dependent)
    end
    vals = promote(A, B, Ḃ, Kₙ, å_default, å)
    T = eltype(vals)
    return HKFActivityModel{T}(vals[1:5]..., vals[6], temperature_dependent)
end

"""
    activity_model(cs::ChemicalSystem, model::HKFActivityModel) -> Function

Return a closure `lna(n, p) -> Vector` computing log-activities for the
extended Debye-Hückel (B-dot) model of Helgeson (1969).

The closure captures all species indices and ionic radii at construction time.
Inside `lna`:
- Solutes: molality convention, B-dot formula for ions, salting-out for neutrals.
- Solvent: osmotic coefficient from Gibbs-Duhem (σ-function).
- Crystals: `ln a = 0` (pure solid).
- Gas: ideal mixture `ln a = ln(xᵢ)`.

If `model.temperature_dependent=true`, `p` must contain `T` (K) and `P` (Pa)
— both are provided automatically by `_build_params`.

AD-compatible: all closure computations accept `ForwardDiff.Dual` inputs.
"""
function activity_model(cs::ChemicalSystem, model::HKFActivityModel)

    # ── Precompute at closure-construction time ────────────────────────────
    idx_solvent = only(cs.idx_solvent)
    idx_solutes = cs.idx_solutes
    idx_gas = cs.idx_gas

    ss_groups = cs.ss_groups
    has_ss = !isempty(ss_groups)
    has_gas = !isempty(idx_gas)
    ss_models = has_ss ? map(ss -> ss.model, cs.solid_solutions) : nothing

    M_w = ustrip(us"kg/mol", cs.species[idx_solvent][:M])   # kg/mol, e.g. 0.018015

    A_fixed = model.A
    B_fixed = model.B
    temp_dep = model.temperature_dependent

    # Per-species data (Float64 — not differentiated).
    zv = Int8[charge(sp) for sp in cs.species]
    åv = Float64[
        iszero(zv[i]) ? 0.0 : _hkf_lookup_å(cs.species[i], model)
            for i in eachindex(zv)
    ]
    n_sp = lastindex(zv)

    idx_ions = [i for i in idx_solutes if !iszero(zv[i])]
    idx_neutrals = [i for i in idx_solutes if  iszero(zv[i])]
    # Setschenow coefficient per neutral species, resolved once (see
    # `_setschenow`): `sp[:Kₙ]` when set, the model's global value otherwise.
    Kₙv = Float64[
        iszero(zv[i]) ? _setschenow(cs.species[i], model) : 0.0 for i in eachindex(zv)
    ]

    ln10 = log(10.0)

    function lna(n::AbstractVector, p)
        ϵ = p.ϵ
        _n = max.(n, ϵ)

        # ── A and B (fixed or T,P-dependent) ──────────────────────────────
        if temp_dep && hasproperty(p, :T) && hasproperty(p, :P)
            AB = hkf_debye_huckel_params(p.T, p.P)
            A, B = AB.A, AB.B
        else
            A, B = A_fixed, B_fixed
        end

        out = zeros(eltype(_n), n_sp)

        # ── Molality: mᵢ = nᵢ / (n_w × M_w) [mol/kg] ─────────────────────
        n_w = _n[idx_solvent]
        denom_mol = n_w * M_w             # kg of solvent

        # ── Ionic strength I = ½ Σ mⱼ zⱼ² ────────────────────────────────
        I = zero(eltype(_n))
        @inbounds for i in idx_solutes
            mᵢ = _n[i] / denom_mol
            I = I + mᵢ * zv[i]^2
        end
        I = I / 2
        sqrtI = sqrt(I + ϵ)              # regularized to avoid Dual NaN at I=0

        # ── Effective radius å_eff for osmotic coefficient ─────────────────
        sum_mz2a = zero(eltype(_n))
        @inbounds for i in idx_ions
            sum_mz2a = sum_mz2a + (_n[i] / denom_mol) * zv[i]^2 * åv[i]
        end
        sum_mz2 = 2 * I                 # Σ mⱼ zⱼ² = 2I by definition
        # Smooth blend: avoids branching on Dual values at ionic-strength ≈ 0.
        # The ϵ term only sets the value in the I → 0 limit, where the osmotic
        # coefficient is 1 regardless; it uses the imposed radius when there is
        # one so that `å = 0` really means a vanishing `B å √I` everywhere.
        å_fallback = model.å === nothing ? model.å_default : model.å
        å_eff = (sum_mz2a + å_fallback * ϵ) / (sum_mz2 + ϵ)

        # ── Ion log-activity coefficients ──────────────────────────────────
        # The formula lives in `_log10γ_ion`, which `activity_coefficients`
        # also calls, so the accessor cannot drift from the solver.
        @inbounds for i in idx_ions
            log10γᵢ = _log10γ_ion(model, zv[i], åv[i], I, sqrtI, A, B)
            mᵢ = _n[i] / denom_mol
            out[i] = ln10 * log10γᵢ + log(mᵢ + ϵ)
        end

        # ── Neutral solute log-activities ──────────────────────────────────
        @inbounds for i in idx_neutrals
            log10γᵢ = _log10γ_neutral(model, I, Kₙv[i])
            mᵢ = _n[i] / denom_mol
            out[i] = ln10 * log10γᵢ + log(mᵢ + ϵ)
        end

        # ── Water activity via osmotic coefficient (Gibbs-Duhem) ───────────
        sum_m = zero(eltype(_n))
        @inbounds for i in idx_solutes
            sum_m = sum_m + _n[i] / denom_mol
        end
        x_arg = B * å_eff * sqrtI
        σ = _hkf_sigma(x_arg)
        φ = 1 - (A * ln10 / 3) * (sum_mz2 / (sum_m + ϵ)) * sqrtI * σ +
            (model.Ḃ * ln10 / 2) * I
        out[idx_solvent] = -M_w * sum_m * φ

        # ── Gas: ideal mixture ─────────────────────────────────────────────
        if has_gas
            n_gas = sum((_n[i] for i in idx_gas); init = zero(eltype(_n)))
            @inbounds for i in idx_gas
                out[i] = log(_n[i] / n_gas)
            end
        end

        # ── Solid solutions ────────────────────────────────────────────────
        if has_ss
            T_val = hasproperty(p, :T) ? p.T : 298.15
            _solid_solution_lna!(out, _n, ss_groups, ss_models, T_val, ϵ)
        end

        return out
    end

    return lna
end

# ── DaviesActivityModel ───────────────────────────────────────────────────────

"""
    struct DaviesActivityModel{T<:Real} <: AbstractActivityModel

The Davies equation ([Davies1962](@cite)): Debye-Hückel with the ion size taken
out, so that nothing per-species has to be known.

# Formulas

```math
\\log_{10}\\gamma_i = -A z_i^2\\left(\\frac{\\sqrt{I}}{1+\\sqrt{I}} - b\\,I\\right) ,
\\qquad
\\log_{10}\\gamma_i = b_n I \\quad (z_i = 0)
```

`\\ln a_i = \\ln 10 \\cdot \\log_{10}\\gamma_i + \\ln m_i` on the molality scale, and
for the solvent **Raoult's law**, `\\ln a_w = \\ln x_w`.

# The physics, and the one thing to know before using it

Davies replaces `1 + B å √I` with `1 + √I`, which amounts to fixing one ion size
for everything (about 3 Å at 25 °C), and adds `b I` to bend the curve back up.
The gain is that **no per-species datum is needed**, which is why it is the model
to reach for when a species list contains ions no radius table covers. The loss
is that every ion of the same charge is now identical, so the model cannot
distinguish Na⁺ from K⁺ at all.

!!! warning "Its water activity does not come from its own activity coefficients"
    The solutes are non-ideal and the solvent is treated as ideal:
    `a_w = x_w` counts molecules. So this model does **not** satisfy the
    Gibbs-Duhem relation between its own `γᵢ` and its `a_w` — the two halves are
    not derived from one excess Gibbs energy. In a dilute solution the error is
    negligible, because both are near their ideal values anyway. In a hydrating
    cement paste, where the water activity is what decides how far the reaction
    goes, the two halves disagree and the answer inherits the disagreement.

    Use [`HKFActivityModel`](@ref), whose `a_w` comes from an osmotic
    coefficient, whenever the water activity is part of the question.

# Inputs, their defaults, and where each default comes from

| field | default | unit | provenance |
|:--|:--|:--|:--|
| `A` | 0.5114 | (kg/mol)^½ | [Helgeson1981](@cite) Table 1 at 25 °C / 1 bar, and reproduced by [`hkf_debye_huckel_params`](@ref) from this package's water model |
| `b` | 0.3 | kg/mol | **part of the published equation** — Davies fixed it, it is not a free parameter of this implementation |
| `bₙ` | 0.1 | kg/mol | a generic salting-out coefficient for neutral species. **No source recorded in this package** |
| `temperature_dependent` | `false` | — | recompute `A` from `p.T`, `p.P` at every call |

# Valid range

`I ≲ 0.5 mol/kg` for the activity coefficients — Davies is usually quoted as
useful to 0.1 and tolerable to 0.5 — and `I ≲ 0.1 mol/kg` for anything that
depends on the water activity.

# Examples
```julia
state_eq = equilibrate(state; model = DaviesActivityModel())

# Temperature-dependent A
state_eq = equilibrate(state; model = DaviesActivityModel(temperature_dependent = true))
```

See also: [`HKFActivityModel`](@ref), [`DiluteSolutionModel`](@ref),
[Activity models](@ref sec-activity-models).
"""
struct DaviesActivityModel{T <: Real} <: AbstractActivityModel
    A::T
    b::T
    bₙ::T
    temperature_dependent::Bool
end

"""
    DaviesActivityModel(; A=0.5114, b=0.3, bₙ=0.1, temperature_dependent=false)

Construct a [`DaviesActivityModel`](@ref).
"""
function DaviesActivityModel(;
        A::Real = 0.5114,
        b::Real = 0.3,
        bₙ::Real = 0.1,
        temperature_dependent::Bool = false,
    )
    vals = promote(A, b, bₙ)
    return DaviesActivityModel{eltype(vals)}(vals..., temperature_dependent)
end

# ── Concentration scale and activity-coefficient formulas ────────────────────
#
# Each model is asked two things beyond its log-activity closure: which
# concentration scale its solute standard state uses, and what its activity
# coefficient is. Both are needed by `activity_coefficients`, `pH(state, model)`
# and the rest of `aqueous_properties.jl`, and the second is called from inside
# the log-activity closures as well, so the formula exists in exactly one place
# and the accessor cannot drift from the solver.

"""
    concentration_scale(model::AbstractActivityModel) -> Symbol

The concentration scale of the model's solute standard state, `:molality` or
`:molarity`.

An activity coefficient is only defined relative to a scale: `a = γ m` on the
molality scale, `a = γ c / c°` on the molarity scale. Nothing in the numbers
says which one a given model used — [`DiluteSolutionModel`](@ref) is on the
molarity scale but takes ρ = 1 kg/L, so its activities coincide numerically with
molalities — so the scale has to be asked rather than inferred. It stops being a
formality as soon as the density departs from 1 kg/L, and it is the origin of
the 0.0013 pH offset that model carries.

See also: [`activity_coefficients`](@ref), [`molalities`](@ref).
"""
concentration_scale(::DiluteSolutionModel) = :molarity
concentration_scale(::HKFActivityModel) = :molality
concentration_scale(::DaviesActivityModel) = :molality

# `log10γ_ion(model, z, å, I, sqrtI, A, B)` and `log10γ_neutral(model, I)` are
# the models' activity-coefficient formulas, as scalars. `å`, `A` and `B` are
# passed in even where a model ignores them, so that all three share one
# signature. AD-safe: arithmetic only.

@inline function _log10γ_ion(model::HKFActivityModel, z, å, I, sqrtI, A, B)
    return -A * z^2 * sqrtI / (1 + B * å * sqrtI) + model.Ḃ * I
end
@inline _log10γ_neutral(model::HKFActivityModel, I) = model.Kₙ * I
# Per-species variant: `Kₙᵢ I` with the species' own coefficient.
@inline _log10γ_neutral(model::HKFActivityModel, I, Kₙᵢ) = Kₙᵢ * I

@inline function _log10γ_ion(model::DaviesActivityModel, z, å, I, sqrtI, A, B)
    return -A * z^2 * (sqrtI / (1 + sqrtI) - model.b * I)
end
@inline _log10γ_neutral(model::DaviesActivityModel, I) = model.bₙ * I

# The dilute model is ideal *on its own scale*: γ = 1 by construction. It is the
# conversion from molality to molarity that makes its activities differ from a
# molality-scale ideal model, not an activity coefficient.
@inline _log10γ_ion(::DiluteSolutionModel, z, å, I, sqrtI, A, B) = zero(I * A)
@inline _log10γ_neutral(::DiluteSolutionModel, I) = zero(I)


"""
    activity_model(cs::ChemicalSystem, model::DaviesActivityModel) -> Function

Return a closure `lna(n, p) -> Vector` computing log-activities for the
Davies (1962) model. No species-specific ionic radii are required.

AD-compatible: all closure computations accept `ForwardDiff.Dual` inputs.
"""
function activity_model(cs::ChemicalSystem, model::DaviesActivityModel)

    idx_solvent = only(cs.idx_solvent)
    idx_solutes = cs.idx_solutes
    idx_gas = cs.idx_gas

    ss_groups = cs.ss_groups
    has_ss = !isempty(ss_groups)
    has_gas = !isempty(idx_gas)
    ss_models = has_ss ? map(ss -> ss.model, cs.solid_solutions) : nothing

    M_w = ustrip(us"kg/mol", cs.species[idx_solvent][:M])

    A_fixed = model.A
    temp_dep = model.temperature_dependent

    zv = Int8[charge(sp) for sp in cs.species]
    n_sp = lastindex(zv)
    idx_ions = [i for i in idx_solutes if !iszero(zv[i])]
    idx_neutrals = [i for i in idx_solutes if  iszero(zv[i])]

    ln10 = log(10.0)

    function lna(n::AbstractVector, p)
        ϵ = p.ϵ
        _n = max.(n, ϵ)

        A = if temp_dep && hasproperty(p, :T) && hasproperty(p, :P)
            hkf_debye_huckel_params(p.T, p.P).A
        else
            A_fixed
        end

        out = zeros(eltype(_n), n_sp)

        n_w = _n[idx_solvent]
        denom_mol = n_w * M_w

        # Ionic strength
        I = zero(eltype(_n))
        @inbounds for i in idx_solutes
            I = I + (_n[i] / denom_mol) * zv[i]^2
        end
        I = I / 2
        sqrtI = sqrt(I + ϵ)

        # Ions — `_log10γ_ion` is shared with `activity_coefficients`.
        @inbounds for i in idx_ions
            log10γᵢ = _log10γ_ion(model, zv[i], 0.0, I, sqrtI, A, 0.0)
            mᵢ = _n[i] / denom_mol
            out[i] = ln10 * log10γᵢ + log(mᵢ + ϵ)
        end

        # Neutral solutes
        @inbounds for i in idx_neutrals
            mᵢ = _n[i] / denom_mol
            out[i] = ln10 * _log10γ_neutral(model, I) + log(mᵢ + ϵ)
        end

        # Water activity — Raoult (mole fraction) approximation
        # n_aqueous ≥ ϵ > 0 always (because _n[i] ≥ ϵ), so no iszero guard needed
        n_aqueous = n_w + sum((_n[i] for i in idx_solutes); init = zero(eltype(_n)))
        out[idx_solvent] = log(n_w / n_aqueous)

        # Gas: ideal mixture
        if has_gas
            n_gas = sum((_n[i] for i in idx_gas); init = zero(eltype(_n)))
            @inbounds for i in idx_gas
                out[i] = log(_n[i] / n_gas)
            end
        end

        # Solid solutions
        if has_ss
            T_val = hasproperty(p, :T) ? p.T : 298.15
            _solid_solution_lna!(out, _n, ss_groups, ss_models, T_val, ϵ)
        end

        return out
    end

    return lna
end

# ── Solid solution activity helpers ───────────────────────────────────────────

"""
    _excess_ln_gamma(model, k, x, T) -> Real

Return the excess log-activity coefficient `ln γₖ` for end-member `k` (1-based index)
of a solid solution with mole-fraction vector `x` at temperature `T` (K).

AD-compatible: all branches preserve `ForwardDiff.Dual` through computations on `x`.

Methods:
- [`IdealSolidSolutionModel`](@ref): returns `zero(eltype(x))`.
- [`RedlichKisterModel`](@ref): binary Redlich-Kister formula (requires `length(x) == 2`).
"""
_excess_ln_gamma(::IdealSolidSolutionModel, k::Int, x::AbstractVector, T::Real) =
    zero(eltype(x))

# Symmetric multi-component Margules. `ln γ_k` is the partial molar derivative
# of `n G^ex / RT` with `G^ex = Σ_{i<j} W_ij x_i x_j`, which gives
# `Σ_{j≠k} W_kj x_j − Σ_{i<j} W_ij x_i x_j`. For two end-members this reduces to
# `W₁₂ x₂²` and `W₁₂ x₁²`, i.e. `RedlichKisterModel(a0 = W₁₂)`.
function _excess_ln_gamma(m::RegularSolutionModel, k::Int, x::AbstractVector, T::Real)
    RT = 8.31446261815324 * T   # J/mol
    n = length(x)
    W = m.W
    lin = zero(eltype(x))
    @inbounds for j in 1:n
        j == k && continue
        lin = lin + (W[k, j] / RT) * x[j]
    end
    quad = zero(eltype(x))
    @inbounds for i in 1:n, j in (i + 1):n
        quad = quad + (W[i, j] / RT) * x[i] * x[j]
    end
    return lin - quad
end

function _excess_ln_gamma(m::RedlichKisterModel, k::Int, x::AbstractVector, T::Real)
    x1, x2 = x[1], x[2]
    RT = 8.31446261815324 * T   # J/mol
    a0 = m.a0 / RT
    a1 = m.a1 / RT
    a2 = m.a2 / RT
    if k == 1
        return x2^2 * (a0 + a1 * (3 * x1 - x2) + a2 * (x1 - x2) * (5 * x1 - x2))
    else
        return x1^2 * (a0 - a1 * (3 * x2 - x1) + a2 * (x2 - x1) * (5 * x2 - x1))
    end
end

"""
    _solid_solution_lna!(out, _n, ss_groups, ss_models, T, ϵ)

Fill `out[i]` with `ln aᵢ = ln xᵢ + ln γᵢ` for all solid-solution end-members.

`ss_groups[k]` and `ss_models[k]` describe the k-th solid-solution phase.
`T` is the temperature in K (only relevant for non-ideal models).
`ϵ` is a regularization floor to avoid `log(0)`.

ForwardDiff-compatible.
"""
function _solid_solution_lna!(
        out::AbstractVector, _n::AbstractVector{ET},
        ss_groups::Vector{Vector{Int}}, ss_models, T, ϵ
    ) where {ET}
    for (grp, mdl) in zip(ss_groups, ss_models)
        n_total = sum(_n[i] for i in grp) + ϵ
        # ET follows eltype(_n) — AD-compatible (Dual when differentiating)
        x = Vector{ET}(undef, length(grp))
        @inbounds for (j, i) in enumerate(grp)
            x[j] = _n[i] / n_total
        end
        @inbounds for (k, i) in enumerate(grp)
            out[i] = log(x[k] + ϵ) + _excess_ln_gamma(mdl, k, x, T)
        end
    end
    return out
end
