# Water retention: the relation between how much water is left in the pore
# space and how tightly it is held.
#
# This is the constitutive input the capillary coupling needs, and it is not
# chemistry. It is measured — by a desorption isotherm, or derived from a pore
# size distribution — it is a property of the particular material, and it evolves
# as that material hydrates. So nothing here ships with a default value. A named
# law whose parameters are not supplied cannot be constructed, and the error
# comes from Julia itself before any number is computed.

using DynamicQuantities

"""
    abstract type WaterRetention

A relation between the degree of saturation of the pore space and the activity
of the water left in it.

The interface is one method: a subtype is **callable**, `r(S) -> a_w`, mapping a
degree of saturation `S ∈ [0, 1]` to a water activity in `(0, 1]`. `S = 1` is a
saturated pore space, where the water is held by nothing and `a_w = 1`.

Two ways to obtain one:

  - a **measured desorption isotherm**, through [`TabulatedRetention`](@ref);
  - a **fitted form** — [`VanGenuchten`](@ref) — whose parameters the caller
    supplies.

And one way to build your own: a pore size distribution plus
[`kelvin_activity`](@ref).

!!! warning "No defaults, deliberately"
    Every named law takes its parameters as keyword arguments **without default
    values**, so omitting one raises `UndefKeywordError` on the spot. A retention
    curve is a measurement on a specific material at a specific age; a default
    would be a number invented on the user's behalf and then silently believed.

See also: [`CapillaryWater`](@ref), [`kelvin_activity`](@ref).
"""
abstract type WaterRetention end

"""
    kelvin_activity(r; γ, V_m, T) -> Real

Activity of water held in a pore of radius `r` by its own meniscus,
`a_w = exp(-2 γ V_m / (r R T))`.

This is the whole of the physics the capillary coupling adds. Water in a fine
pore is at a lower chemical potential than bulk water of the same composition,
so hydrates that consume water become less stable, and below some pore size the
reaction stops — which is what self-desiccation is.

Every parameter is a keyword **without a default**: the surface tension `γ`
belongs to the liquid and the temperature, and the molar volume `V_m` to the
liquid, so neither is the package's to assume.

# Examples

```jldoctest
julia> using DynamicQuantities

julia> a = kelvin_activity(4.2e-9u"m"; γ = 0.072u"N/m", V_m = 1.8e-5u"m^3/mol", T = 298.15u"K");

julia> round(a; digits = 2)      # the self-desiccation plateau of a sealed paste
0.78
```

See also: [`kelvin_radius`](@ref), [`WaterRetention`](@ref).
"""
function kelvin_activity(r; γ, V_m, T)
    r_m = ustrip(us"m", r)
    r_m > 0 || throw(ArgumentError("kelvin_activity: the pore radius must be positive, got $r"))
    return exp(
        -2 * ustrip(us"N/m", γ) * ustrip(us"m^3/mol", V_m) /
            (r_m * ustrip(us"J/mol/K", Constants.R) * ustrip(us"K", T))
    )
end

"""
    kelvin_radius(a_w; γ, V_m, T) -> Real

The pore radius, in meters, whose meniscus holds water at activity `a_w`. The
inverse of [`kelvin_activity`](@ref).

Useful for reading a retention curve back as a pore size, which is the form in
which the literature usually discusses it: `a_w = 0.78` — the relative humidity a
sealed high-performance paste settles at — is 4.2 nm, the gel-pore scale.

# Examples

```jldoctest
julia> using DynamicQuantities

julia> r = kelvin_radius(0.78; γ = 0.072u"N/m", V_m = 1.8e-5u"m^3/mol", T = 298.15u"K");

julia> round(r * 1e9; digits = 1)      # nanometres
4.2
```

See also: [`kelvin_activity`](@ref).
"""
function kelvin_radius(a_w; γ, V_m, T)
    0 < a_w < 1 || throw(
        ArgumentError(
            "kelvin_radius: the water activity must lie strictly between 0 and 1, " *
                "got $a_w — at a_w = 1 the meniscus is flat and the radius infinite"
        )
    )
    return -2 * ustrip(us"N/m", γ) * ustrip(us"m^3/mol", V_m) /
        (log(a_w) * ustrip(us"J/mol/K", Constants.R) * ustrip(us"K", T))
end

"""
    TabulatedRetention(; S, a_w) -> TabulatedRetention

A measured desorption isotherm, as two vectors: degrees of saturation and the
water activities (relative humidities) at which they were observed.

Interpolation is linear in `ln a_w` against `S`, which is the variable the
capillary coupling actually uses and the one in which a Kelvin curve is closest
to straight. Both keywords are mandatory.

`S` must be strictly increasing and `a_w` non-decreasing with it — a retention
curve that falls as the pore space fills is not a retention curve, and an inner
constructor refuses it rather than interpolating nonsense. Outside the tabulated
range the value is clamped to the nearest end and a warning is issued once: an
isotherm measured between `S = 0.3` and `S = 1` says nothing about `S = 0.05`,
and extrapolating a logarithm there produces confident absurdity.

# Examples

```julia
# RH measured at four saturations on a hardened paste
r = TabulatedRetention(; S = [0.30, 0.50, 0.75, 1.00], a_w = [0.44, 0.66, 0.85, 1.00])
r(0.60)      # interpolated
```

See also: [`WaterRetention`](@ref), [`VanGenuchten`](@ref).
"""
struct TabulatedRetention{T <: Real} <: WaterRetention
    S::Vector{T}
    a_w::Vector{T}
    warned::Base.RefValue{Bool}

    function TabulatedRetention{T}(S::AbstractVector, a_w::AbstractVector) where {T <: Real}
        length(S) == length(a_w) || throw(
            ArgumentError(
                "TabulatedRetention: S and a_w must have the same length, got " *
                    "$(length(S)) and $(length(a_w))."
            )
        )
        length(S) >= 2 || throw(
            ArgumentError("TabulatedRetention: at least two points are needed to interpolate.")
        )
        all(>(0), diff(S)) || throw(
            ArgumentError("TabulatedRetention: S must be strictly increasing.")
        )
        all(>=(0), diff(a_w)) || throw(
            ArgumentError(
                "TabulatedRetention: a_w must not decrease as S increases — a pore " *
                    "space that holds its water more tightly as it fills is not a " *
                    "retention curve."
            )
        )
        all(x -> 0 < x <= 1, a_w) || throw(
            ArgumentError("TabulatedRetention: every a_w must lie in (0, 1].")
        )
        return new{T}(Vector{T}(S), Vector{T}(a_w), Ref(false))
    end
end

TabulatedRetention(; S, a_w) =
    TabulatedRetention{float(promote_type(eltype(S), eltype(a_w)))}(S, a_w)

function (r::TabulatedRetention)(S::Real)
    x = r.S
    if S < x[1] || S > x[end]
        if !r.warned[]
            @warn "TabulatedRetention: the saturation left the tabulated range; " *
                "the value is clamped to the nearest measured end. An isotherm says " *
                "nothing outside the saturations it was measured at." S range = (x[1], x[end])
            r.warned[] = true
        end
        return S < x[1] ? r.a_w[1] : r.a_w[end]
    end
    k = searchsortedlast(x, S)
    k >= length(x) && return r.a_w[end]
    t = (S - x[k]) / (x[k + 1] - x[k])
    # Linear in ln a_w: the coupling reads the logarithm, and a Kelvin curve is
    # closest to straight there.
    return exp((1 - t) * log(r.a_w[k]) + t * log(r.a_w[k + 1]))
end

"""
    VanGenuchten(; a, m) -> VanGenuchten

The van Genuchten water-retention form, as a capillary pressure

```
p_c(S) = a * (S^(-1/m) - 1)^(1 - m)
```

converted to a water activity by the Kelvin relation
`a_w = exp(-p_c V_m / (R T))`, with `V_m` and `T` supplied by the constraint
that uses it.

Both parameters are mandatory keywords: `a` is a capillary pressure scale (Pa,
or a pressure quantity) and `m` is the dimensionless shape exponent, `0 < m < 1`.
There are no defaults — the pair belongs to a fitted material, not to the
package.

# Published parameters, and the two conventions

[BaroghelBouny1999](@cite) fit exactly this expression to measured water-vapor
desorption isotherms, and write it with `b = 1/m`:

```
p_c(S) = a (S^(-b) - 1)^(1 - 1/b)
```

so a `b` from that literature becomes `m = 1/b` here. Their Table 5, for
materials whose mixes and porosities are in their Tables 1 and 4:

| mix | material | W/C | `a` (MPa) | `b` | `m = 1/b` |
|:--|:--|:--|:--|:--|:--|
| CO | cement paste | 0.34 | 37.5479 | 2.1684 | 0.46117 |
| CH | paste, 10 % silica fume | 0.19 | 96.2837 | 1.9540 | 0.51177 |
| BO | concrete | 0.48 | 18.6237 | 2.2748 | 0.43960 |
| BH | concrete, 10 % silica fume | 0.26 | 46.9364 | 2.0601 | 0.48541 |

# Examples

```julia
# The ordinary cement paste of Baroghel-Bouny et al. (1999), their mix CO
r = VanGenuchten(; a = 37.5479e6, m = 1 / 2.1684)
```

!!! warning "A pair belongs to a material, not to this package"
    The table above is quoted with its source so it can be checked, not so it
    can be copied blindly. Those four fits are one cement, one curing history and
    one age; `m` and `a` both move with w/c, with silica fume and with the
    aggregate. Fit your own isotherm where you have one, and cite what you used.

See also: [`WaterRetention`](@ref), [`TabulatedRetention`](@ref).
"""
struct VanGenuchten{T <: Real} <: WaterRetention
    a::T
    m::T

    function VanGenuchten{T}(a, m) where {T <: Real}
        a > 0 || throw(ArgumentError("VanGenuchten: the pressure scale `a` must be positive, got $a."))
        0 < m < 1 || throw(
            ArgumentError("VanGenuchten: the exponent `m` must lie strictly between 0 and 1, got $m.")
        )
        return new{T}(a, m)
    end
end

function VanGenuchten(; a, m)
    a_pa = a isa Number && !(a isa DynamicQuantities.AbstractQuantity) ? float(a) :
        ustrip(us"Pa", a)
    vals = promote(a_pa, float(m))
    return VanGenuchten{eltype(vals)}(vals...)
end

"""
    capillary_pressure(r::VanGenuchten, S) -> Real

The capillary pressure in pascals at saturation `S`, before the Kelvin relation
turns it into an activity. Exposed because it is the form the retention
literature tabulates and fits.
"""
function capillary_pressure(r::VanGenuchten, S::Real)
    Sc = clamp(S, 1.0e-12, 1.0)
    Sc >= 1 && return zero(r.a)
    return r.a * (Sc^(-1 / r.m) - 1)^(1 - r.m)
end

# A `VanGenuchten` is not callable on `S` alone: turning its pressure into an
# activity needs the liquid's molar volume and the temperature, which the
# constraint holds. `_retention_activity` is the internal hook the constraint
# calls, with a generic fallback for the laws that ARE a plain function of `S`.
_retention_activity(r::WaterRetention, S, V_m, T) = r(S)
function _retention_activity(r::VanGenuchten, S, V_m, T)
    return exp(-capillary_pressure(r, S) * V_m / (ustrip(us"J/mol/K", Constants.R) * T))
end

# A bare function is a retention law too — the escape hatch for a caller whose
# curve is neither tabulated nor van Genuchten.
struct FunctionRetention{F} <: WaterRetention
    f::F
end
(r::FunctionRetention)(S::Real) = r.f(S)
