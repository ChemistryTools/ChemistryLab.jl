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
    solvent_fraction(state) -> Float64

Mole fraction of the aqueous solvent within the aqueous phase,
`n_w / Σ_aqueous n`. One for pure water, and the number that says whether there
is still a solution to speak of.

Every quantity built on the aqueous phase — molality, ionic strength, activity,
pH — is defined *per kilogram of solvent*, and the dual solver parameterizes its
interior variables by the solvent's chemical potential. All of that presumes the
solvent is the phase, not one species in it. A real electrolyte, even a
concentrated one, keeps `x_w` above about 0.9: seawater is 0.99, a saturated NaCl
brine 0.90.

Below [`SOLVENT_FRACTION_FLOOR`](@ref) the formulation has no ground left.
Measured on a sealed cement paste driven under its stoichiometric water demand —
w/c = 0.28 on the mix of the [w/c example](@ref sec-wc-ratio) — the Gibbs
minimum consumes the free water down to 6e-9 mol, `x_w` falls to 0.21, and the
ionic strength comes out at 409 mol/kg against a Debye-Huckel model valid to
about one. The amounts such a solve reports are not equilibrium values; the
answer is that the question was posed outside the model's domain.

See also: [`molalities`](@ref), [`ionic_strength`](@ref),
[`equilibrate_certified`](@ref), which checks this on the answer it returns.
"""
function solvent_fraction(state::ChemicalState)
    cs = state.system
    isempty(cs.idx_solvent) && return 0.0
    i_w = only(cs.idx_solvent)
    n = ustrip.(us"mol", state.n)
    tot = sum(_primal(n[i]) for i in cs.idx_aqueous; init = 0.0)
    return tot <= 0 ? 0.0 : _primal(n[i_w]) / tot
end

"""
    SOLVENT_FRACTION_FLOOR

The mole fraction of solvent below which an aqueous phase is no longer a
solution, and every quantity derived from it is meaningless. See
[`solvent_fraction`](@ref).

`0.5` is not a modeling choice but a floor no real solution comes near: it is
about 28 mol of solute per kilogram of water, past saturation for anything. A
value below it means the solve has driven the water into the solids, not that the
solution is concentrated.
"""
const SOLVENT_FRACTION_FLOOR = 0.5

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
    ionic_strength(state::ChemicalState; kind = :effective, ϵ = 1e-16) -> Float64

Molality-basis ionic strength `I = ½ Σⱼ mⱼ zⱼ²`, in mol/kg.

`kind` selects which of the two ionic strengths is meant, and they are different
quantities:

  - `:effective` (the default) sums over the **speciated free ions** — the
    composition as solved, with every complex and ion pair counted at its own
    charge. A neutral pair such as `Ca(SO4)@` contributes nothing. This is the
    one every activity model in this package is a function of, and the one to
    pass to any of them.
  - `:stoichiometric` sums as if every complex were **fully dissociated** into
    the components of `state.system`, so `Ca(SO4)@` contributes as
    `Ca²⁺ + SO₄²⁻`. This is the analytical ionic strength of the recipe rather
    than of the solution, and it is what some correlations for salting-out and
    for diffusivity are fitted against.

The gap between them measures how much of the salt is associated. On a Portland
cement pore solution it is small — 0.2121 against 0.2136 mol/kg, 0.7 % — because
little is paired at that ionic strength; on a sulfate brine it is not.

Neither is a property of the activity model: both are properties of the
composition. Compare the effective one before comparing anything else — an
ionic strength that disagrees between two codes means they are not describing
the same solution, whatever their volumes happen to agree on.

The stoichiometric sum uses the decomposition of each species over the system's
primaries (`state.system.SM.A`), counting `|νₚ| zₚ²` for each charged primary
`p`. That is the same decomposition the mass balance uses, so it needs no
separate table of dissociation reactions; note that it attributes `OH⁻`, which
CEMDATA18 writes as `H₂O − H⁺`, one unit of charge through the `H⁺` component,
which is the right count.

Throws if the system has no aqueous phase, or on an unknown `kind`.

# Examples

```julia
ionic_strength(eq)                          # 0.2121 mol/kg — free ions
ionic_strength(eq; kind = :stoichiometric)  # 0.2136 — fully dissociated
```

See also: [`molalities`](@ref), [`activity_coefficients`](@ref).
"""
function ionic_strength(
        state::ChemicalState; kind::Symbol = :effective, ϵ::Float64 = 1.0e-16
    )
    cs = state.system
    _require_aqueous(cs, "ionic_strength")
    m = molalities(state; ϵ = ϵ)

    if kind === :effective
        I = 0.0
        for i in cs.idx_solutes
            z = Int(charge(cs.species[i]))
            iszero(z) && continue
            I += m[symbol(cs.species[i])] * z^2
        end
        return I / 2
    elseif kind === :stoichiometric
        A = cs.SM.A
        prim = cs.SM.primaries
        # `Zz` is the charge pseudo-component of CEMDATA18, not an ion.
        charged = [
            (p, Int(charge(prim[p]))^2) for p in eachindex(prim)
                if symbol(prim[p]) != "Zz" && !iszero(charge(prim[p]))
        ]
        I = 0.0
        for i in cs.idx_solutes
            mi = m[symbol(cs.species[i])]
            iszero(mi) && continue
            for (p, z2) in charged
                I += mi * abs(Float64(A[p, i])) * z2
            end
        end
        return I / 2
    end
    throw(
        ArgumentError(
            "ionic_strength: kind must be :effective or :stoichiometric, got :$kind",
        )
    )
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
        log10γ = if iszero(z)
            # Same per-species Setschenow coefficient the closure uses.
            model isa HKFActivityModel ?
                _log10γ_neutral(model, I, _setschenow(cs.species[i], model)) :
                _log10γ_neutral(model, I)
        else
            _log10γ_ion(model, z, åv[i], I, sqrtI, AB.A, AB.B)
        end
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

# ── Automatic initial approximation by continuation ──────────────────────────
#
# A Gibbs energy minimization needs a starting point, and on a cement the choice
# decides whether it converges at all. From the "cold" state — all the mass in
# the reactants, every product at the `ϵ` floor — no back end reaches the
# optimum on a CEM I paste of 135 species: the answer comes back
# `optimal = false` with a worst supersaturation of order 1e1 and a total volume
# 13 % wrong. The problem is convex, so this is not a local minimum; it is the
# conditioning of an interior-point method started against the boundary. With
# 120 of 135 species at 1e-16 and nine at ~1 mol, the barrier gradients span
# sixteen orders of magnitude and the fraction-to-the-boundary rule crawls; and
# a solid-solution end-member has `ln a = ln x → −∞` as its mole fraction goes
# to zero, so the objective's gradient is unbounded on exactly the face where a
# mixing phase vanishes.
#
# GEM-Selektor solves this by computing an initial approximation rather than
# asking for one: `AutoInitialApproximation` (GEMS3K, `ipm_simplex.cpp`) is an
# "LPP-based automatic initial approximation of the primal vector x" obtained
# with a "modified simplex method with two-side constraints". Reproducing that
# needs a genuine LP solver — a barrier method on a linear objective is not one,
# and measured here it does not move off the cold state at all.
#
# What does work, and needs no new algorithm, is continuation. Scale everything
# but the solvent by λ and walk λ from a small value to 1, each step started
# from the answer to the previous one. At small λ the reactants are dilute, the
# solid amounts are small, and the problem sits far from the awkward face; the
# assemblage then grows continuously with λ. Measured on the CEM I paste: the
# first rung does not certify and it does not matter — every rung from λ = 0.02
# on certifies, and the endpoint is identical, to the digits printed, to the
# answer from a chemically informed seed.
#
# Nothing here is chemical. "Everything but the solvent" needs no knowledge of
# which hydrates will form, and λ is not a physical parameter: at λ = 1 the
# state is exactly the one given.

"""
    saturation_indices(state, model; ϵ = 1e-16) -> OrderedDict{String, <:Real}

`LogSI` for every species at `state`: `log₁₀(IAP/K)` of the reaction that forms it
from the system's primaries.

Zero for a phase at equilibrium with the solution, negative for one that is
undersaturated, positive for one that **should have precipitated**. It is the
quantity GEM-Selektor prints as `LogSI`, and the one an
[`optimality_certificate`](@ref) summarizes into a single worst violation without
saying which phase that is.

No fitting is involved. The row labels of the conservation matrix are the primary
species, so a component's element potential is that primary's chemical potential,
`y_c = μ_c/RT`, and

```
LogSI_s = [Σ_c A_cs y_c − μ_s/RT] / ln 10
```

which is [`saturation_ratio`](@ref) written for the formation reaction
`s = Σ_c A_cs (primary c)`. A conservation row that labels no species — the
charge row — contributes nothing, since the coefficient of any neutral phase
there is zero.

Two things to know before reading the numbers:

  - **The check is built in.** Every phase actually present at an equilibrium
    must come out at `LogSI = 0`; measured on a CEM I paste, the twelve present
    solids land within 1.2e-12. If they do not, the state is not an equilibrium
    and no other index in the result means anything.
  - **A solid-solution end-member's index is relative to its current mole
    fraction**, since its activity is `ln x`. For an end-member at the solver's
    lower bound that is a statement about a vanishing phase, not about whether
    the solid solution would form.

# Examples

```julia
si = saturation_indices(eq, model)
si["hydrotalcite"]                     # +5.58: absent, and it should not be
[k for (k, v) in si if v > 1e-4]       # everything supersaturated
```

See also: [`optimality_certificate`](@ref), [`saturation_ratio`](@ref),
[`log_activities`](@ref).
"""
function saturation_indices(
        state::ChemicalState, model::AbstractActivityModel = DiluteSolutionModel();
        ϵ::Float64 = 1.0e-16,
    )
    cs = state.system
    lna = log_activities(state, model; ϵ = ϵ)
    p = _build_params(state; ϵ = ϵ)
    g = [p.ΔₐG⁰overT[i] + lna[symbol(cs.species[i])] for i in eachindex(cs.species)]
    A = cs.SM.A
    idx = Dict(symbol(sp) => i for (i, sp) in enumerate(cs.species))
    # A row whose primary is not among the species — the charge row — gets zero,
    # which is exact for every neutral phase.
    y = [get(idx, symbol(pr), 0) for pr in cs.SM.primaries]
    yv = [k == 0 ? zero(eltype(g)) : g[k] for k in y]
    inv_ln10 = inv(log(10))
    return OrderedDict(
        symbol(cs.species[i]) =>
            (sum(A[c, i] * yv[c] for c in eachindex(yv)) - g[i]) * inv_ln10
            for i in eachindex(cs.species)
    )
end

"""
    homotopy_initial_state(state::ChemicalState; model = DiluteSolutionModel(),
                           steps = ..., ϵ = 1e-16, max_bisections = 6,
                           balance_atol = 1e-5, balance_rtol = 1e-3,
                           verbose = false)
        -> Union{ChemicalState, Nothing}

An automatic initial approximation obtained by continuation in the amount of
solute.

Scales every species except the aqueous solvent by a factor `λ`, and walks `λ`
through `steps` up to 1, solving at each value from the answer to the previous
one. At `λ = 1` the composition is exactly `state`, so the returned composition
respects the same element balance; it is a **starting point**, not a certified
equilibrium.

`steps` is a suggestion, not a schedule. A rung is **accepted only if it
conserves mass**, and a refused rung is retaken by halving the distance back to
the last `λ` that worked, up to `max_bisections` times per target. Both halves of
that matter: a rung that has wandered off the balance would otherwise be carried
forward as the start of every later rung, and a rung that simply cannot be taken
in one jump can be taken in two.

The acceptance test is `|rᵢ| ≤ balance_atol + balance_rtol · scaleᵢ` on the
element-balance residual, row by row, with
`scaleᵢ = max(|bᵢ|, Σⱼ |Aᵢⱼ| nⱼ)`. It has to be mixed rather than relative: two
conservation rows carry a legitimately negligible budget — electroneutrality is
exactly zero, and a cement recipe is routinely given a carbon trace of 1e-9 mol
— so a relative test rejects residuals of 1e-10 mol as if they were failures.
Neither tolerance certifies anything; `balance_atol` sits above the accuracy the
interior point itself reaches (~3e-6 mol here) because a rung is a guess, and
the certificate judges the result afterwards.

This is what makes a realistic cement solvable from a cold start: with all the
mass in the reactants and every product at the `ϵ` floor, no back end reaches
the optimum, because an interior-point method started against the boundary of a
system spanning sixteen orders of magnitude in amount cannot take a useful step.
At small `λ` the same system is dilute and well away from that boundary.

It requires **nothing** from the caller beyond the initial state: no guess at
which phases will form, no knowledge of the answer. [`equilibrate_certified`](@ref)
calls it automatically when its ordinary starting points fail to certify, so
ordinary use never needs it.

The walk is done under `model`, which defaults to [`DiluteSolutionModel`](@ref)
and should normally be left there **even when the target is a non-ideal model**.
Measured on the CEM I paste: walking under the extended Debye-Hückel model with
a common ion size of zero runs away to an ionic strength of 18 mol/kg, because
that model's coefficients fall steeply with `I`, which raises solubility, which
raises `I`. The ideal model has no such feedback, walks cleanly, and its answer
is a good starting point for the non-ideal one — which is how
[`equilibrate_certified`](@ref) uses it.

Returns `nothing` if the walk produced nothing usable.

Differentiability is unaffected. Under `ForwardDiff`, `equilibrate_certified`
strips to the primal state, solves in `Float64` and attaches the sensitivity
through the implicit function theorem, so this never sees a `Dual` and the
derivative does not depend on how the starting point was found.

# Examples

```julia
# What `equilibrate_certified` does for you when the cold start fails:
guess = homotopy_initial_state(state)
eq, cert = equilibrate_certified(guess; model = HKFActivityModel())
```

See also: [`equilibrate_certified`](@ref).
"""
function homotopy_initial_state(
        state::ChemicalState;
        model::AbstractActivityModel = DiluteSolutionModel(),
        steps = (0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.85, 1.0),
        ϵ::Float64 = 1.0e-16,
        max_bisections::Int = 6,
        balance_atol::Float64 = 1.0e-5,
        balance_rtol::Float64 = 1.0e-3,
        verbose::Bool = false,
    )
    cs = state.system
    isempty(cs.idx_solvent) && return nothing
    i_w = only(cs.idx_solvent)
    n0 = ustrip.(us"mol", state.n)

    # `STRICT_CONVERGENCE[]` has to be off along the walk, and restored after.
    #
    # The walk *relies* on its early rungs being allowed to fall short: measured
    # on a CEM I paste, the first one (λ = 0.01) does not converge and it does
    # not matter, because all it has to be is close to the next. Under strict
    # convergence that rung raises instead, the step is skipped, every later rung
    # starts cold again, and the continuation degenerates into the very failure
    # it exists to avoid. These are guesses, not answers; the answer is the
    # certified solve that follows, and that one is still judged strictly.
    #
    # `_EXPLORING_STARTS` goes with it, for the same reason one step further: a
    # rung that ends on `MaxIters` is expected, and warning about it makes a walk
    # that worked read as a walk that failed. `verbose = true` reports every rung
    # either way.
    strict = STRICT_CONVERGENCE[]
    STRICT_CONVERGENCE[] = false
    try
        return _exploring_starts() do
            _homotopy_walk(
                cs, i_w, n0, model, steps, ϵ, verbose, max_bisections,
                balance_atol, balance_rtol,
            )
        end
    finally
        STRICT_CONVERGENCE[] = strict
    end
end

"""
    _homotopy_rung(cs, A, i_w, n0, model, λ, start, ϵ, verbose, atol, rtol)
        -> Union{ChemicalState, Nothing}

Solve one rung of the continuation, and **accept it only if it conserves mass**.

`nothing` means "refuse this rung", which the walk answers by taking a smaller
step. Two things can go wrong at a rung, and only one of them raises: a back end
can throw, or it can return an answer that is not on the constraint surface. The
second is the dangerous one, because the walk would then carry that composition
forward as the start of every later rung. Measured on a CEM I paste: an accepted
rung off the balance by moles takes the whole walk with it, and
`equilibrate_certified` ends on an answer with an element balance of 6.7 mol —
every hydrate at zero, and a table of amounts that reads like a result.

The test is the mixed one, `|rᵢ| ≤ atol + rtol · scaleᵢ`, and it has to be mixed.
A purely relative test is meaningless on this problem because two conservation
rows carry a legitimately negligible budget: electroneutrality is **exactly
zero**, and a cement recipe is routinely given a carbon trace of 1e-9 mol.
Measured at λ = 0.01 on that paste, a residual of 1.7e-10 mol on the charge row
scores 168 against its own budget and a residual of 1.1e-10 mol on the carbon
row scores 10.6 — both physically nothing, both rejected, and the walk then
never reached its first three targets while still appearing to work. A row
holding 1e-11 mol cannot be balanced better than the solver's absolute floor,
and asking it to be is a category error.

`scaleᵢ = max(|bλᵢ|, Σⱼ |Aᵢⱼ| nⱼ)`: the row's own budget, or the amount of that
component actually being moved around when the budget is near zero — the natural
yardstick for a conservation row that nets to nothing. `atol` sits above the
accuracy the interior point itself reaches (about 3e-6 mol on this class of
problem), because a rung is a guess and not an answer. Neither number certifies
anything: they are there to reject a rung that has *wandered*, and the
certificate judges the result afterwards.
"""
function _homotopy_rung(cs, A, i_w, n0, model, λ, start, ϵ, verbose, atol, rtol)
    nλ = [i == i_w ? n0[i] : λ * n0[i] for i in eachindex(n0)]
    bλ = A * nλ
    from = start === nothing ? ChemicalState(cs, nλ .* u"mol") : start
    for f in _SOLVER_FACTORIES
        stepped = nothing
        try
            esolver = EquilibriumSolver(cs, model, f())
            stepped = SciMLBase.solve(esolver, from; ϵ = ϵ, b = bλ)
        catch err
            verbose && @info "homotopy rung raised" λ = λ backend = f err
            continue
        end
        n = Float64[ustrip(us"mol", x) for x in stepped.n]
        scale = max.(abs.(bλ), abs.(A) * n)
        # Dimensionless, and 1 is the boundary: the largest residual measured
        # against the tolerance allowed for its own row.
        off = maximum(abs.(A * n - bλ) ./ (atol .+ rtol .* scale))
        if off <= 1
            return stepped
        end
        verbose && @info "homotopy rung off the balance" λ = λ backend = f excess = off
    end
    return nothing
end

function _homotopy_walk(
        cs, i_w, n0, model, steps, ϵ, verbose, max_bisections, atol, rtol,
    )
    A = Float64.(cs.SM.A)
    current = nothing        # the answer at `done`
    done = 0.0               # the largest λ actually reached
    for target in steps
        λ = target
        # Aim for the target; on a refused rung, halve the distance back to the
        # last λ that worked and try again, and after a rung that lands, aim for
        # the target once more from there. This is what makes the walk robust
        # rather than lucky: the fixed ladder is a suggestion, and a rung that
        # cannot be taken in one jump is taken in two. Bounded, so a target that
        # cannot be reached at all costs a handful of solves and is skipped.
        for _ in 0:max_bisections
            done >= target && break
            stepped = _homotopy_rung(
                cs, A, i_w, n0, model, λ, current, ϵ, verbose, atol, rtol,
            )
            if stepped === nothing
                λ = 0.5 * (done + λ)
                verbose && @info "homotopy step halved" λ = λ
            else
                current, done = stepped, λ
                verbose && @info "homotopy step" λ = λ
                λ = target
            end
        end
        done >= target ||
            verbose && @info "homotopy target not reached" target = target reached = done
    end
    return current
end
