# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using TOML

"""
    build_pitzer_parameters(toml_file) -> PitzerParameters

Read a Pitzer interaction-parameter set from a TOML file.

Nothing about this is automatic: the file has to be named, which is the whole
point. `data/pitzer-reardon1990.toml` ships with the package and is reached as

```julia
p = build_pitzer_parameters(datapath("pitzer-reardon1990.toml"))
```

# File format

```toml
[meta]
name = "..."
temperature_K = 298.15
speciation = "dissociated"        # informational; read by the caller, not here

[[binary]]                        # one per cation-anion pair
cation = "Na+"
anion  = "Cl-"
beta0  = 0.0765
beta1  = 0.2664
beta2  = 0.0
Cphi   = 0.0013
origin = "fitted"                 # or "estimated:HSO4-"

[[theta]]                         # like-charge pair, either order
i = "Na+"
j = "K+"
value = -0.012

[[psi]]                           # triplet, any order
i = "Na+"
j = "Cl-"
k = "SO4-2"
value = 0.0014

[[lambda]]                        # neutral with an ion
neutral = "H2CO3@"
ion = "Na+"
value = 0.1
origin = "fitted"
```

`beta1`, `beta2` and `Cphi` default to zero within a `[[binary]]` entry;
`beta0` does not, because an entry without it describes nothing. The shape
constants `alpha1`, `alpha1_22`, `alpha2` and `b` may be given under `[meta]`
and otherwise take the conventional values documented in
[`PitzerParameters`](@ref).

# What it refuses

A `[[binary]]` entry without `cation`, `anion` or `beta0`, and a duplicate pair
— silently keeping the last of two conflicting entries would be worse than
stopping. Species names are **not** checked against any database here: whether
the set covers the system is decided when the model is attached to it, which is
where the error can name the pairs that are missing.

See also: [`PitzerParameters`](@ref), [`PitzerActivityModel`](@ref),
[`pitzer_origin`](@ref).
"""
function build_pitzer_parameters(toml_file::AbstractString)
    path = resolve_data_path(toml_file)
    raw = TOML.parsefile(path)
    meta = get(raw, "meta", Dict{String, Any}())

    beta0 = Dict{Tuple{String, String}, Float64}()
    beta1 = Dict{Tuple{String, String}, Float64}()
    beta2 = Dict{Tuple{String, String}, Float64}()
    Cphi = Dict{Tuple{String, String}, Float64}()
    origin = Dict{Tuple{String, String}, String}()

    for (i, e) in enumerate(get(raw, "binary", Any[]))
        for key in ("cation", "anion", "beta0")
            haskey(e, key) || error(
                "build_pitzer_parameters: [[binary]] entry $i in $path has no \"$key\"."
            )
        end
        k = (String(e["cation"]), String(e["anion"]))
        haskey(beta0, k) && error(
            "build_pitzer_parameters: $path gives the pair $(k[1])/$(k[2]) twice."
        )
        beta0[k] = float(e["beta0"])
        beta1[k] = float(get(e, "beta1", 0.0))
        beta2[k] = float(get(e, "beta2", 0.0))
        Cphi[k] = float(get(e, "Cphi", 0.0))
        origin[k] = String(get(e, "origin", "unrecorded"))
    end

    theta = Dict{Tuple{String, String}, Float64}()
    for (i, e) in enumerate(get(raw, "theta", Any[]))
        haskey(e, "i") && haskey(e, "j") && haskey(e, "value") || error(
            "build_pitzer_parameters: [[theta]] entry $i in $path needs i, j and value."
        )
        theta[(String(e["i"]), String(e["j"]))] = float(e["value"])
    end

    psi = Dict{Tuple{String, String, String}, Float64}()
    for (i, e) in enumerate(get(raw, "psi", Any[]))
        all(haskey(e, k) for k in ("i", "j", "k", "value")) || error(
            "build_pitzer_parameters: [[psi]] entry $i in $path needs i, j, k and value."
        )
        psi[(String(e["i"]), String(e["j"]), String(e["k"]))] = float(e["value"])
    end

    lambda = Dict{Tuple{String, String}, Float64}()
    for (i, e) in enumerate(get(raw, "lambda", Any[]))
        all(haskey(e, k) for k in ("neutral", "ion", "value")) || error(
            "build_pitzer_parameters: [[lambda]] entry $i in $path needs " *
                "neutral, ion and value."
        )
        k = (String(e["neutral"]), String(e["ion"]))
        lambda[k] = float(e["value"])
        haskey(e, "origin") && (origin[k] = String(e["origin"]))
    end

    return PitzerParameters(;
        beta0 = beta0, beta1 = beta1, beta2 = beta2, Cphi = Cphi,
        theta = theta, psi = psi, lambda = lambda,
        alpha1 = float(get(meta, "alpha1", 2.0)),
        alpha1_22 = float(get(meta, "alpha1_22", 1.4)),
        alpha2 = float(get(meta, "alpha2", 12.0)),
        b = float(get(meta, "b", 1.2)),
        origin = origin,
    )
end
