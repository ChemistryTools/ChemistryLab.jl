# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# The Pitzer ion-interaction model. Unlike the extended Debye-Hückel models in
# `activities.jl`, this one is not a per-species formula: the excess Gibbs energy
# is a virial expansion, so a coefficient belongs to a *pair* of ions and another
# to a *triplet*, and the sums run over the whole solution. That is why it needs
# its own closure rather than a method of `_log10γ_ion`, and why its parameters
# are caller input: no thermodynamic database ships them.

using DynamicQuantities

"""
    PitzerParameters(; beta0, beta1, beta2, Cphi, theta, psi, lambda,
                       alpha1_22 = 1.4, alpha1 = 2.0, alpha2 = 12.0, b = 1.2,
                       origin = Dict{Tuple{String, String}, String}())

The interaction parameters of a Pitzer model, keyed by species symbol.

Every keyword is **mandatory**, `origin` and the four shape constants excepted:
a Pitzer model is only as good as the set it is given, and a default would be a
number invented on the caller's behalf. Omitting one raises `UndefKeywordError`
before anything is computed.

# The tables

| keyword | keys | meaning |
|:--|:--|:--|
| `beta0`, `beta1`, `beta2` | `(cation, anion)` | the second virial coefficient of that pair, as three terms of one ionic-strength dependence |
| `Cphi` | `(cation, anion)` | the third virial coefficient of that pair |
| `theta` | `(ion, ion)` of **like** charge | interaction between two cations, or two anions |
| `psi` | `(ion, ion, ion)` | a triplet: two anions with a cation, or two cations with an anion |
| `lambda` | `(neutral, ion)` | a neutral solute with an ion — salting out, in Pitzer's form |

`theta`, `psi` and `lambda` entries are looked up in any order of their keys.
A **missing** `theta`, `psi` or `lambda` is taken as zero, which is the
convention of the literature these tables come from; a **missing** `beta0` for a
pair the system actually contains is an error, raised when the model is attached
to a [`ChemicalSystem`](@ref) — see [`PitzerActivityModel`](@ref).

# The shape constants

`alpha1 = 2.0`, `alpha1_22 = 1.4` (used when both ions carry a charge of
magnitude 2 or more), `alpha2 = 12.0` and `b = 1.2` kg^½·mol^−½ are **not**
fitted parameters: they fix the functional form, and the published `β` tables
were fitted assuming them. Changing one invalidates the table it is used with.
They are keywords so that a set fitted with a different convention can be used,
not so that they can be tuned.

# Provenance

`origin` maps a pair to a free-text note, and the loader
[`build_pitzer_parameters`](@ref) fills it from the data file. It exists because
a published set can contain values that were *estimated by analogy* rather than
measured, and a caller reading `beta0` has otherwise no way of telling the two
apart. [`pitzer_origin`](@ref) reads it back.

See also: [`PitzerActivityModel`](@ref), [`build_pitzer_parameters`](@ref).
"""
struct PitzerParameters{T <: Real}
    beta0::Dict{Tuple{String, String}, T}
    beta1::Dict{Tuple{String, String}, T}
    beta2::Dict{Tuple{String, String}, T}
    Cphi::Dict{Tuple{String, String}, T}
    theta::Dict{Tuple{String, String}, T}
    psi::Dict{Tuple{String, String, String}, T}
    lambda::Dict{Tuple{String, String}, T}
    alpha1::T
    alpha1_22::T
    alpha2::T
    b::T
    origin::Dict{Tuple{String, String}, String}

    # Inner, so that the generated constructor cannot bypass the checks: a
    # `beta1` for a pair with no `beta0` is a typo in the caller's table, and a
    # `theta` between ions of opposite charge is a category error — that
    # interaction is what `beta` is for.
    function PitzerParameters{T}(
            beta0, beta1, beta2, Cphi, theta, psi, lambda,
            alpha1, alpha1_22, alpha2, b, origin,
        ) where {T <: Real}
        for (name, tbl) in (("beta1", beta1), ("beta2", beta2), ("Cphi", Cphi))
            for k in keys(tbl)
                haskey(beta0, k) || throw(
                    ArgumentError(
                        "PitzerParameters: $name has an entry for $k but beta0 does not. " *
                            "A pair is described by all four coefficients or by none.",
                    )
                )
            end
        end
        b > 0 || throw(ArgumentError("PitzerParameters: b must be positive, got $b."))
        alpha1 > 0 && alpha2 > 0 && alpha1_22 > 0 || throw(
            ArgumentError("PitzerParameters: the alpha constants must be positive.")
        )
        return new{T}(
            beta0, beta1, beta2, Cphi, theta, psi, lambda,
            alpha1, alpha1_22, alpha2, b, origin,
        )
    end
end

function PitzerParameters(;
        beta0, beta1, beta2, Cphi, theta, psi, lambda,
        alpha1::Real = 2.0, alpha1_22::Real = 1.4, alpha2::Real = 12.0, b::Real = 1.2,
        origin::AbstractDict = Dict{Tuple{String, String}, String}(),
    )
    T = float(
        promote_type(
            eltype(values(beta0)), eltype(values(beta1)), eltype(values(beta2)),
            eltype(values(Cphi)), eltype(values(theta)), eltype(values(psi)),
            eltype(values(lambda)), typeof(alpha1), typeof(alpha1_22),
            typeof(alpha2), typeof(b),
        )
    )
    conv(d, K) = Dict{K, T}(k => T(v) for (k, v) in d)
    P2 = Tuple{String, String}
    P3 = Tuple{String, String, String}
    return PitzerParameters{T}(
        conv(beta0, P2), conv(beta1, P2), conv(beta2, P2), conv(Cphi, P2),
        conv(theta, P2), conv(psi, P3), conv(lambda, P2),
        T(alpha1), T(alpha1_22), T(alpha2), T(b),
        Dict{P2, String}(k => String(v) for (k, v) in origin),
    )
end

"""
    pitzer_origin(p::PitzerParameters, cation, anion) -> String

What the parameters of that pair are: `"fitted"`, `"estimated:<analog>"`, or
`"unrecorded"` when the set carries no note.

A published Pitzer set can contain values obtained by analogy with a chemically
similar ion rather than by fitting a measurement — Reardon's cement set does so
for every silicate, aluminate and ferrate pair, which are exactly the ions a
cement assemblage needs. This accessor is how that shows up in a calculation
instead of staying in a paper's footnote.
"""
pitzer_origin(p::PitzerParameters, cation::AbstractString, anion::AbstractString) =
    get(p.origin, (String(cation), String(anion)), "unrecorded")

# Symmetric lookups: the caller's table may key a pair or a triplet either way
# round, and the physics does not distinguish them.
@inline function _sym2(d, i, j)
    v = get(d, (i, j), nothing)
    v === nothing || return v
    return get(d, (j, i), zero(valtype(d)))
end
@inline function _sym3(d, i, j, k)
    for key in ((i, j, k), (i, k, j), (j, i, k), (j, k, i), (k, i, j), (k, j, i))
        v = get(d, key, nothing)
        v === nothing || return v
    end
    return zero(valtype(d))
end

"""
    PitzerActivityModel(; parameters, temperature_dependent = false)

The Pitzer ion-interaction activity model.

`parameters` is a [`PitzerParameters`](@ref) and has **no default**: this model
cannot be reached without one, which is the point. Whether the set is complete
is not a property of the set alone but of the set *relative to a species list*,
so the check happens when the model meets a [`ChemicalSystem`](@ref): every
cation-anion pair the system contains must have a `beta0` entry, and the error
names the pairs that do not.

# Why this model rather than an extended Debye-Hückel one

`γ` and the osmotic coefficient are derived from **one** excess Gibbs energy — a
virial expansion in the molalities, with one coefficient per ion pair and one
per triplet — so the Gibbs-Duhem relation between the solutes and the solvent
holds by construction rather than approximately. That is what
[`DaviesActivityModel`](@ref) does not do and [`HKFActivityModel`](@ref) does
only up to a mean-radius approximation. See
[Activity models](@ref sec-theory-activity).

# What is *not* implemented

  - the **higher-order electrostatic terms** ``{}^E\\theta(I)`` and
    ``{}^E\\theta'(I)``, which correct the interaction of two ions of
    *unsymmetrical* charge (a 1+ with a 2+, say). They are exactly zero for a
    symmetrical pair, so a single 1-1 or 2-2 electrolyte is unaffected; in a
    mixture of Na⁺ and Ca²⁺ they are a real omission, and the docstring says so
    rather than the code pretending otherwise;
  - any temperature dependence of the interaction parameters themselves. A
    published set is fitted at one temperature — 25 °C for the set shipped here
    — and this model does not extrapolate it. `temperature_dependent` governs
    only ``A_\\varphi``, which comes from the water model.

# Example

```julia
p = build_pitzer_parameters(datapath("pitzer-reardon1990.toml"))
model = PitzerActivityModel(; parameters = p)
```

See also: [`PitzerParameters`](@ref), [`build_pitzer_parameters`](@ref),
[`pitzer_origin`](@ref).
"""
struct PitzerActivityModel{T <: Real} <: AbstractActivityModel
    parameters::PitzerParameters{T}
    temperature_dependent::Bool
end

PitzerActivityModel(; parameters, temperature_dependent::Bool = false) =
    PitzerActivityModel(parameters, temperature_dependent)

concentration_scale(::PitzerActivityModel) = :molality

# ── The two ionic-strength functions of the model ────────────────────────────
#
# `g` and `g′` are 0/0 at zero ionic strength, so both carry a Taylor branch.
# The series are exact to the order written: with e^{-x} expanded,
# `1 - (1+x)e^{-x} = x²/2 - x³/3 + x⁴/8`, hence `g = 1 - 2x/3 + x²/4`, and
# `1 - (1+x+x²/2)e^{-x} = x³/6 - x⁴/8`, hence `g′ = -x/3 + x²/4`.
# Branching is on the primal so that a `ForwardDiff.Dual` keeps its derivative.
@inline function _pitzer_g(x::T) where {T <: Real}
    if abs(_primal(x)) < 1.0e-4
        return one(T) - (2 // 3) * x + x^2 / 4
    end
    return 2 * (one(T) - (one(T) + x) * exp(-x)) / x^2
end

@inline function _pitzer_gp(x::T) where {T <: Real}
    if abs(_primal(x)) < 1.0e-4
        return -x / 3 + x^2 / 4
    end
    return -2 * (one(T) - (one(T) + x + x^2 / 2) * exp(-x)) / x^2
end

"""
    activity_model(cs::ChemicalSystem, model::PitzerActivityModel) -> Function

Return the closure `lna(n, p)` of the Pitzer model for `cs`.

The completeness of `model.parameters` is checked **here**, against the species
`cs` actually contains, and a missing cation-anion pair raises rather than
defaulting to ideal behavior.
"""
function activity_model(cs::ChemicalSystem, model::PitzerActivityModel)
    idx_solvent = only(cs.idx_solvent)
    idx_solutes = cs.idx_solutes
    idx_gas = cs.idx_gas
    ss_groups = cs.ss_groups
    has_ss = !isempty(ss_groups)
    has_gas = !isempty(idx_gas)
    ss_models = has_ss ? map(ss -> ss.model, cs.solid_solutions) : nothing

    M_w = ustrip(us"kg/mol", cs.species[idx_solvent][:M])
    n_sp = lastindex(cs.species)
    par = model.parameters

    sym = [symbol(sp) for sp in cs.species]
    zv = Int8[charge(sp) for sp in cs.species]
    cats = [i for i in idx_solutes if zv[i] > 0]
    ans = [i for i in idx_solutes if zv[i] < 0]
    neus = [i for i in idx_solutes if iszero(zv[i])]

    # ── The refusal that makes "no half-parameterized Pitzer" true ───────────
    missing_pairs = Tuple{String, String}[]
    for c in cats, a in ans
        haskey(par.beta0, (sym[c], sym[a])) || push!(missing_pairs, (sym[c], sym[a]))
    end
    if !isempty(missing_pairs)
        listed = join(("$(c)/$(a)" for (c, a) in missing_pairs), ", ")
        throw(
            ArgumentError(
                "PitzerActivityModel: the parameter set has no beta0 for " *
                    "$(length(missing_pairs)) cation-anion pair(s) this system contains: " *
                    "$listed. A Pitzer model cannot fall back on ideal behavior for a " *
                    "pair it does not describe, so the model refuses rather than " *
                    "returning a number. Either supply the missing parameters or build " *
                    "a species list the set covers. Note that a set fitted for a " *
                    "dissociated speciation does not describe ion pairs such as " *
                    "Ca(SO4)@ or CaOH+: those associations are already inside its beta " *
                    "coefficients, and carrying them as species counts them twice.",
            )
        )
    end

    # Precomputed per-pair tables, Float64 and not differentiated.
    npair = (length(cats), length(ans))
    B0 = [par.beta0[(sym[c], sym[a])] for c in cats, a in ans]
    B1 = [get(par.beta1, (sym[c], sym[a]), 0.0) for c in cats, a in ans]
    B2 = [get(par.beta2, (sym[c], sym[a]), 0.0) for c in cats, a in ans]
    CC = [
        get(par.Cphi, (sym[c], sym[a]), 0.0) /
            (2 * sqrt(abs(Int(zv[c]) * Int(zv[a])))) for c in cats, a in ans
    ]
    A1 = [
        (abs(zv[c]) >= 2 && abs(zv[a]) >= 2) ? par.alpha1_22 : par.alpha1
            for c in cats, a in ans
    ]
    ΘCC = [_sym2(par.theta, sym[i], sym[j]) for i in cats, j in cats]
    ΘAA = [_sym2(par.theta, sym[i], sym[j]) for i in ans, j in ans]
    ΨCCA = [_sym3(par.psi, sym[i], sym[j], sym[a]) for i in cats, j in cats, a in ans]
    ΨAAC = [_sym3(par.psi, sym[i], sym[j], sym[c]) for i in ans, j in ans, c in cats]
    ΛNC = [_sym2(par.lambda, sym[nn], sym[i]) for nn in neus, i in cats]
    ΛNA = [_sym2(par.lambda, sym[nn], sym[i]) for nn in neus, i in ans]

    α2 = par.alpha2
    bp = par.b
    A_fixed = 0.5114                     # log10-basis Debye-Hückel A at 25 °C
    temp_dep = model.temperature_dependent

    function lna(n::AbstractVector, p)
        ϵ = p.ϵ
        _n = max.(n, ϵ)
        TT = eltype(_n)
        out = zeros(TT, n_sp)

        A_log10 = if temp_dep && hasproperty(p, :T) && hasproperty(p, :P)
            hkf_debye_huckel_params(p.T, p.P).A
        else
            A_fixed
        end
        Aφ = A_log10 * log(10) / 3       # the osmotic-coefficient basis

        kg = _n[idx_solvent] * M_w
        m = [_n[i] / kg for i in eachindex(_n)]

        I = zero(TT)
        Σm = zero(TT)
        Z = zero(TT)
        @inbounds for i in idx_solutes
            I += m[i] * zv[i]^2
            Σm += m[i]
            Z += m[i] * abs(zv[i])
        end
        I /= 2
        sI = sqrt(I + ϵ)

        # f^γ and the pair sums that make up F
        fγ = -Aφ * (sI / (1 + bp * sI) + 2 / bp * log1p(bp * sI))
        F = fγ
        @inbounds for (ic, c) in enumerate(cats), (ia, a) in enumerate(ans)
            Bp = (
                B1[ic, ia] * _pitzer_gp(A1[ic, ia] * sI) +
                    B2[ic, ia] * _pitzer_gp(α2 * sI)
            ) / (I + ϵ)
            F += m[c] * m[a] * Bp
        end

        # Σ_c Σ_a m_c m_a C_ca, shared by every ion
        ΣmmC = zero(TT)
        @inbounds for (ic, c) in enumerate(cats), (ia, a) in enumerate(ans)
            ΣmmC += m[c] * m[a] * CC[ic, ia]
        end

        # ── cations ──────────────────────────────────────────────────────────
        @inbounds for (ic, c) in enumerate(cats)
            s = zv[c]^2 * F + abs(zv[c]) * ΣmmC
            for (ia, a) in enumerate(ans)
                B = B0[ic, ia] + B1[ic, ia] * _pitzer_g(A1[ic, ia] * sI) +
                    B2[ic, ia] * _pitzer_g(α2 * sI)
                s += m[a] * (2 * B + Z * CC[ic, ia])
            end
            for (jc, c2) in enumerate(cats)
                c2 == c && continue
                acc = 2 * ΘCC[ic, jc]
                for (ia, a) in enumerate(ans)
                    acc += m[a] * ΨCCA[ic, jc, ia]
                end
                s += m[c2] * acc
            end
            for (ia, a) in enumerate(ans), (ja, a2) in enumerate(ans)
                ja > ia || continue
                s += m[a] * m[a2] * ΨAAC[ia, ja, ic]
            end
            for (inn, nn) in enumerate(neus)
                s += 2 * m[nn] * ΛNC[inn, ic]
            end
            out[c] = s + log(m[c] + ϵ)
        end

        # ── anions ───────────────────────────────────────────────────────────
        @inbounds for (ia, a) in enumerate(ans)
            s = zv[a]^2 * F + abs(zv[a]) * ΣmmC
            for (ic, c) in enumerate(cats)
                B = B0[ic, ia] + B1[ic, ia] * _pitzer_g(A1[ic, ia] * sI) +
                    B2[ic, ia] * _pitzer_g(α2 * sI)
                s += m[c] * (2 * B + Z * CC[ic, ia])
            end
            for (ja, a2) in enumerate(ans)
                a2 == a && continue
                acc = 2 * ΘAA[ia, ja]
                for (ic, c) in enumerate(cats)
                    acc += m[c] * ΨAAC[ia, ja, ic]
                end
                s += m[a2] * acc
            end
            for (ic, c) in enumerate(cats), (jc, c2) in enumerate(cats)
                jc > ic || continue
                s += m[c] * m[c2] * ΨCCA[ic, jc, ia]
            end
            for (inn, nn) in enumerate(neus)
                s += 2 * m[nn] * ΛNA[inn, ia]
            end
            out[a] = s + log(m[a] + ϵ)
        end

        # ── neutral solutes ──────────────────────────────────────────────────
        @inbounds for (inn, nn) in enumerate(neus)
            s = zero(TT)
            for (ic, c) in enumerate(cats)
                s += 2 * m[c] * ΛNC[inn, ic]
            end
            for (ia, a) in enumerate(ans)
                s += 2 * m[a] * ΛNA[inn, ia]
            end
            out[nn] = s + log(m[nn] + ϵ)
        end

        # ── the solvent, from the osmotic coefficient ────────────────────────
        acc = -Aφ * I^(3 // 2) / (1 + bp * sI)
        @inbounds for (ic, c) in enumerate(cats), (ia, a) in enumerate(ans)
            Bφ = B0[ic, ia] + B1[ic, ia] * exp(-A1[ic, ia] * sI) +
                B2[ic, ia] * exp(-α2 * sI)
            acc += m[c] * m[a] * (Bφ + Z * CC[ic, ia])
        end
        @inbounds for (ic, c) in enumerate(cats), (jc, c2) in enumerate(cats)
            jc > ic || continue
            t = ΘCC[ic, jc]
            for (ia, a) in enumerate(ans)
                t += m[a] * ΨCCA[ic, jc, ia]
            end
            acc += m[c] * m[c2] * t
        end
        @inbounds for (ia, a) in enumerate(ans), (ja, a2) in enumerate(ans)
            ja > ia || continue
            t = ΘAA[ia, ja]
            for (ic, c) in enumerate(cats)
                t += m[c] * ΨAAC[ia, ja, ic]
            end
            acc += m[a] * m[a2] * t
        end
        @inbounds for (inn, nn) in enumerate(neus)
            for (ic, c) in enumerate(cats)
                acc += m[nn] * m[c] * ΛNC[inn, ic]
            end
            for (ia, a) in enumerate(ans)
                acc += m[nn] * m[a] * ΛNA[inn, ia]
            end
        end
        φ = 1 + 2 * acc / (Σm + ϵ)
        out[idx_solvent] = -M_w * Σm * φ

        if has_gas
            n_gas = sum((_n[i] for i in idx_gas); init = zero(TT))
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
