# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities
using OrderedCollections

# ── Aqueous properties of a solved state ─────────────────────────────────────
#
# The activity closures compute the molalities, the ionic strength and the
# activity coefficients on their way to the log-activities, and used to keep all
# three to themselves. Anyone comparing a state against GEM-Selektor, PHREEQC or
# Reaktoro needs them species by species, so they are read back here.
#
# Two things this file is careful about, because both have produced wrong
# comparisons:
#
#   * γ is computed from the model's formula, not as a ratio a/m. A species
#     parked at the solver's lower bound has its log-activity dominated by the
#     closures' `+ ϵ` regularization, and the ratio then returns nonsense —
#     values of order 1e300 for a charge class whose only members are trace.
#     The formula is a function of the ionic strength and the charge alone and
#     is exact whatever the amount.
#   * An activity is a number on a scale, and nothing in the number says which.
#     `DiluteSolutionModel` puts its solutes on the molarity scale and takes
#     ρ = 1 kg/L, so its activities happen to equal the molalities as numbers;
#     the other two are on the molality scale, where that equality holds by
#     construction. [`concentration_scale`](@ref) is the only way to tell them
#     apart, and it stops being a formality as soon as ρ departs from 1 kg/L.

# Index of the aqueous solvent, or an error naming the caller. Every function in
# this file needs a solvent to divide by, so failing here is better than
# returning a silently meaningless number.
function _require_aqueous(cs::ChemicalSystem, what::AbstractString)
    isempty(cs.idx_solvent) && throw(
        ArgumentError(
            "$what needs an aqueous phase: no species of class `SC_AQSOLVENT` " *
                "is present in this system.",
        )
    )
    return only(cs.idx_solvent)
end

# Debye-Hückel A and B in the model's own convention. Davies has no B (its
# denominator is 1 + √I), and the dilute model has neither.
_debye_huckel_AB(model::HKFActivityModel, T_K, P_Pa) =
    model.temperature_dependent ? hkf_debye_huckel_params(T_K, P_Pa) :
    (A = model.A, B = model.B)
function _debye_huckel_AB(model::DaviesActivityModel, T_K, P_Pa)
    A = model.temperature_dependent ? hkf_debye_huckel_params(T_K, P_Pa).A : model.A
    return (A = A, B = zero(A))
end
_debye_huckel_AB(::DiluteSolutionModel, T_K, P_Pa) = (A = 0.0, B = 0.0)

# The effective radius actually used for each species, by the same lookup the
# closure uses. Zero for models that have no radius.
function _ion_sizes(cs::ChemicalSystem, model::HKFActivityModel)
    return Float64[
        iszero(charge(sp)) ? 0.0 : _hkf_lookup_å(sp, model) for sp in cs.species
    ]
end
_ion_sizes(cs::ChemicalSystem, ::AbstractActivityModel) = zeros(Float64, length(cs.species))

"""
    molalities(state::ChemicalState; ϵ = 1e-16) -> OrderedDict{String,Float64}

Molality `mᵢ = nᵢ / (n_w Mw)` of every aqueous solute, in mol per kg of solvent.

The solvent itself is not a solute and is omitted. `ϵ` floors the amounts the
same way the activity closures do, so the values match what the solver saw;
species at the floor come back at a molality of order `ϵ` rather than zero.

Throws if the system has no aqueous phase.

# Examples

```julia
m = molalities(eq)
m["K+"]                      # mol/kg of water
sum(values(m))               # total solute molality
```

See also: [`ionic_strength`](@ref), [`activity_coefficients`](@ref),
[`concentration_scale`](@ref).
"""
function molalities(state::ChemicalState; ϵ::Float64 = 1.0e-16)
    cs = state.system
    i_w = _require_aqueous(cs, "molalities")
    n = ustrip.(us"mol", state.n)
    M_w = ustrip(us"kg/mol", cs.species[i_w][:M])
    kg_solvent = max(_primal(n[i_w]), ϵ) * M_w
    out = OrderedDict{String, Float64}()
    for i in cs.idx_solutes
        out[symbol(cs.species[i])] = max(_primal(n[i]), ϵ) / kg_solvent
    end
    return out
end

"""
    ionic_strength(state::ChemicalState; ϵ = 1e-16) -> Float64

Molality-basis ionic strength `I = ½ Σⱼ mⱼ zⱼ²`, in mol/kg.

This is a property of the composition, not of the activity model: every model in
the package computes it this way, and it is the quantity their coefficients are
functions of. Compare it before comparing anything else — an ionic strength that
disagrees means the two codes are not describing the same solution, whatever
their volumes happen to agree on.

Throws if the system has no aqueous phase.

# Examples

```julia
ionic_strength(eq)           # e.g. 0.212 mol/kg for a CEM I pore solution
```

See also: [`molalities`](@ref), [`activity_coefficients`](@ref).
"""
function ionic_strength(state::ChemicalState; ϵ::Float64 = 1.0e-16)
    cs = state.system
    _require_aqueous(cs, "ionic_strength")
    m = molalities(state; ϵ = ϵ)
    I = 0.0
    for i in cs.idx_solutes
        z = Int(charge(cs.species[i]))
        iszero(z) && continue
        I += m[symbol(cs.species[i])] * z^2
    end
    return I / 2
end

"""
    log_activities(state::ChemicalState, model::AbstractActivityModel)
        -> OrderedDict{String,Float64}

Natural log of the activity of **every** species, in the model's own convention.

This is the vector the Gibbs energy is built from: `μᵢ/RT = ΔₐG⁰ᵢ/RT + ln aᵢ`.
Crystals of a pure phase get `ln a = 0`, solid-solution end-members `ln a = ln xᵢ`
plus their excess term, the solvent its osmotic contribution, and solutes the
log of their concentration in the model's scale plus `ln γᵢ`.

A species at the solver's lower bound has its value dominated by the closures'
`+ ϵ` regularization; read [`activity_coefficients`](@ref) rather than dividing
these by a concentration.

# Examples

```julia
lna = log_activities(eq, model)
lna["H2O@"]                                  # ln a_w
exp(lna["Portlandite"])                      # 1.0 for a pure phase that is present
```

See also: [`activities`](@ref), [`activity_coefficients`](@ref).
"""
function log_activities(
        state::ChemicalState, model::AbstractActivityModel; ϵ::Float64 = 1.0e-16
    )
    cs = state.system
    lna_fun = activity_model(cs, model)
    p = _build_params(state; ϵ = ϵ)
    n = ustrip.(us"mol", state.n)
    lna = lna_fun(n, p)
    out = OrderedDict{String, Float64}()
    for (i, sp) in enumerate(cs.species)
        out[symbol(sp)] = _primal(lna[i])
    end
    return out
end

"""
    activities(state::ChemicalState, model::AbstractActivityModel)
        -> OrderedDict{String,Float64}

Activity of every species — `exp` of [`log_activities`](@ref).

# Examples

```julia
a = activities(eq, HKFActivityModel())
a["H2O@"]                    # water activity, e.g. 0.9937
```
"""
function activities(
        state::ChemicalState, model::AbstractActivityModel; ϵ::Float64 = 1.0e-16
    )
    lna = log_activities(state, model; ϵ = ϵ)
    return OrderedDict{String, Float64}(k => exp(v) for (k, v) in lna)
end

"""
    activity_coefficients(state::ChemicalState, model::AbstractActivityModel)
        -> OrderedDict{String,Float64}

Activity coefficient γᵢ of every aqueous species, from the model's formula.

Solutes are evaluated as `γᵢ = 10^(log₁₀ γᵢ)` with the model's own expression —
`−A zᵢ² √I/(1 + B åᵢ √I) + Ḃ I` for the ions of [`HKFActivityModel`](@ref), `Kₙ I`
for its neutrals, and identically 1 for [`DiluteSolutionModel`](@ref), which is
ideal *on its own* (molarity) scale. The solvent is reported as `γ_w = a_w / x_w`.

**Not** computed as a ratio of activity to concentration. That ratio agrees for
an abundant solute — and the tests check that it does — but it diverges for a
species parked at the solver's lower bound, whose log-activity is dominated by
the closures' `+ ϵ` term: it returns values of order 1e300 for a charge class
whose only members are trace. The formula depends on the ionic strength and the
charge alone, so it is exact for a trace species and for a major one alike.

Throws if the system has no aqueous phase.

# Examples

```julia
γ = activity_coefficients(eq, HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0))
γ["K+"], γ["Ca+2"]           # 0.6098, 0.1199 on a CEM I pore solution
γ["H2O@"]                    # the solvent, as a_w / x_w
```

See also: [`ionic_strength`](@ref), [`concentration_scale`](@ref),
[`activities`](@ref).
"""
function activity_coefficients(
        state::ChemicalState, model::AbstractActivityModel; ϵ::Float64 = 1.0e-16
    )
    cs = state.system
    i_w = _require_aqueous(cs, "activity_coefficients")
    n = ustrip.(us"mol", state.n)

    I = ionic_strength(state; ϵ = ϵ)
    sqrtI = sqrt(I)
    T_K = ustrip(us"K", temperature(state))
    P_Pa = ustrip(us"Pa", pressure(state))
    AB = _debye_huckel_AB(model, T_K, P_Pa)
    åv = _ion_sizes(cs, model)

    out = OrderedDict{String, Float64}()
    for i in cs.idx_solutes
        z = Int(charge(cs.species[i]))
        log10γ = iszero(z) ? _log10γ_neutral(model, I) :
            _log10γ_ion(model, z, åv[i], I, sqrtI, AB.A, AB.B)
        out[symbol(cs.species[i])] = 10.0^_primal(log10γ)
    end

    # The solvent has no formula of that shape: its activity comes from the
    # osmotic coefficient (HKF) or from Raoult (the other two), so report the
    # coefficient that the mole-fraction convention implies.
    n_aq = sum(max(_primal(n[i]), ϵ) for i in cs.idx_aqueous)
    x_w = max(_primal(n[i_w]), ϵ) / n_aq
    a_w = exp(_primal(log_activities(state, model; ϵ = ϵ)[symbol(cs.species[i_w])]))
    out[symbol(cs.species[i_w])] = a_w / x_w

    return out
end

"""
    pH(state::ChemicalState, model::AbstractActivityModel) -> Float64

`−log₁₀ a(H⁺)` — the pH in the **activity** convention of `model`.

This is what GEM-Selektor and Reaktoro report, and it is **not** what the
one-argument [`pH`](@ref) returns: that one is `−log₁₀ c(H⁺)` with the
concentration taken over the computed liquid volume, and in an alkaline solution
it is reconstructed from OH⁻ through `pKw`. The two differ by the activity
coefficient and by the scale conversion. On a Portland cement pore solution at
I ≈ 0.2 mol/kg, with γ(H⁺) ≈ 0.61, the gap is about **0.21 units** — large
enough to be mistaken for a modeling error when comparing against another code.

Returns `NaN` if the system carries no `H+`.

# Examples

```julia
pH(eq)                       # 13.310 — concentration convention
pH(eq, model)                # 13.099 — activity convention, comparable to GEMS
```

See also: [`pOH`](@ref), [`activity_coefficients`](@ref), [`FixedpH`](@ref).
"""
function pH(
        state::ChemicalState, model::AbstractActivityModel; ϵ::Float64 = 1.0e-16
    )
    return _p_activity(state, model, "H+"; ϵ = ϵ)
end

"""
    pOH(state::ChemicalState, model::AbstractActivityModel) -> Float64

`−log₁₀ a(OH⁻)` — the pOH in the **activity** convention of `model`.

See [`pH`](@ref) for why this differs from the one-argument [`pOH`](@ref).

Returns `NaN` if the system carries no `OH-`.
"""
function pOH(
        state::ChemicalState, model::AbstractActivityModel; ϵ::Float64 = 1.0e-16
    )
    return _p_activity(state, model, "OH-"; ϵ = ϵ)
end

function _p_activity(
        state::ChemicalState, model::AbstractActivityModel, sym::AbstractString;
        ϵ::Float64 = 1.0e-16,
    )
    cs = state.system
    _require_aqueous(cs, "pH(state, model)")
    i = findfirst(s -> symbol(s) == sym, cs.species)
    i === nothing && return NaN
    lna = log_activities(state, model; ϵ = ϵ)
    return -lna[sym] / log(10)
end
