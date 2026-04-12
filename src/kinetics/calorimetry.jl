using DynamicQuantities
using LinearAlgebra

# ── AbstractCalorimeter ───────────────────────────────────────────────────────

"""
    abstract type AbstractCalorimeter end

Base type for calorimeter models that can be coupled to a kinetics simulation.

A calorimeter augments the ODE state vector `u` with additional variables
(cumulative heat Q and/or temperature T) and extends the ODE right-hand side
accordingly.

Concrete subtypes:
- [`IsothermalCalorimeter`](@ref): T = constant, tracks Q(t) = ∫q̇dt.
- [`SemiAdiabaticCalorimeter`](@ref): dT/dt = (q̇ - L(T-T_env)) / Cp.

All methods must be AD-compatible (no `Float64` casts).
"""
abstract type AbstractCalorimeter end

# ── Heat-rate from kinetic reactions ─────────────────────────────────────────

"""
    heat_rate(kinetic_reactions, rates, T_K; ϵ=1e-30) -> Real

Compute the instantaneous heat generation rate [W = J/s]:

```
q̇ = Σᵢ rᵢ(t) × ΔHᵣ,ᵢ(T)
```

where `rᵢ` [mol/s] is the net rate of the i-th kinetic reaction (positive =
dissolution/forward) and `ΔHᵣ,ᵢ(T)` [J/mol] is the enthalpy of reaction
at temperature `T_K`.

`ΔHᵣ` is computed as the sum of stoichiometrically weighted species enthalpies
`ΔₐH⁰` using the thermodynamic functions already stored in each species.
If a species lacks the `:ΔₐH⁰` property the contribution is treated as zero.

# Arguments

  - `kinetic_reactions`: vector of [`KineticReaction`](@ref).
  - `rates`: vector of net rates [mol/s] (length = number of kinetic reactions).
  - `T_K`: temperature [K] (plain number or Quantity; stripped internally).
  - `ϵ`: regularisation floor (unused here, kept for API uniformity).

AD-compatible: the ΔₐH⁰ callables accept `ForwardDiff.Dual` T inputs.
"""
function heat_rate(
        kinetic_reactions::AbstractVector,
        rates::AbstractVector,
        T_K;
        kwargs...,
    )
    T_val = ustrip(T_K)
    q = zero(promote_type(eltype(rates), typeof(T_val)))

    for (kr, r) in zip(kinetic_reactions, rates)
        ΔHr = _reaction_enthalpy(kr, T_val)
        q = q + r * ΔHr
    end
    return q
end

# ── _reaction_enthalpy dispatch hierarchy ────────────────────────────────────
#
# Priority (highest to lowest):
#   1. KineticReaction with explicit heat_per_mol (Float64) — user-supplied value
#   2. KineticReaction with heat_per_mol=nothing — delegate to reaction stoichiometry
#   3. AbstractReaction — stoichiometric sum of species ΔₐH⁰
#   4. AbstractSpecies — no stoichiometry available, returns 0

# 1. Explicit heat_per_mol: user-supplied enthalpy of reaction [J/mol].
#    Bypasses stoichiometry — correct for empirical models (e.g. ParrotKilloh).
function _reaction_enthalpy(kr::KineticReaction{<:Any, <:Any, <:Any, Float64}, ::Real)
    return kr.heat_per_mol
end

# 2. No explicit heat_per_mol: delegate to the stored reaction object.
function _reaction_enthalpy(kr::KineticReaction{<:Any, <:Any, <:Any, Nothing}, T_K::Real)
    return _reaction_enthalpy(kr.reaction, T_K)
end

# 3. Full reaction stoichiometry: ΔHᵣ = Σ_products ν_j·ΔₐH⁰_j - Σ_reactants ν_i·ΔₐH⁰_i
function _reaction_enthalpy(reaction::AbstractReaction, T_K::Real)
    ΔH = zero(T_K)
    for (sp, ν) in pairs(reactants(reaction))
        haskey(properties(sp), :ΔₐH⁰) || continue
        ΔH -= ν * ustrip(sp[:ΔₐH⁰](T = T_K * u"K"; unit = true))
    end
    for (sp, ν) in pairs(products(reaction))
        haskey(properties(sp), :ΔₐH⁰) || continue
        ΔH += ν * ustrip(sp[:ΔₐH⁰](T = T_K * u"K"; unit = true))
    end
    return ΔH
end

# 4. Bare species (convenience constructor, no reaction stoichiometry): return 0.
#    For meaningful calorimetry with bare-species KineticReactions, supply
#    `heat_per_mol` to the KineticReaction constructor.
_reaction_enthalpy(::AbstractSpecies, T_K::Real) = zero(T_K)

# ── Total-enthalpy helper ─────────────────────────────────────────────────────

"""
    _total_enthalpy(n_full, h_fns, T_K) -> Real

Compute the total molar enthalpy of the system:

```math
H = \\sum_i n_i \\cdot \\Delta_a H^\\circ_i(T)
```

where `n_full` is the full species mole vector [mol] and `h_fns` is the vector of
standard-enthalpy callables built by [`build_kinetics_params`](@ref) (one entry per
species, `nothing` for species without `:ΔₐH⁰` data).

Called by the `DiscreteCallback` in `KineticsOrdinaryDiffEqExt` at each accepted ODE
step to track the total system enthalpy. The cumulative heat released is then:

```math
Q(t) = H(0) - H(t)
```

This captures heat from **all** chemical transformations — both kinetically controlled
reactions and instantaneous equilibrium re-speciation.

AD-compatible: accepts `ForwardDiff.Dual` numbers in `n_full` and `T_K`.
"""
function _total_enthalpy(n_full::AbstractVector, h_fns, T_K::Real)
    H = zero(promote_type(eltype(n_full), typeof(T_K)))
    for (i, hf) in enumerate(h_fns)
        isnothing(hf) && continue
        h_i = ustrip(hf(; T = T_K * u"K", unit = true))
        H += n_full[i] * h_i
    end
    return H
end

# ── IsothermalCalorimeter ─────────────────────────────────────────────────────

"""
    struct IsothermalCalorimeter{T<:Real} <: AbstractCalorimeter

Isothermal calorimeter: temperature is held constant at `T` [K], and the
cumulative heat released `Q(t) = ∫₀ᵗ q̇(τ) dτ` [J] is tracked as an
additional ODE state.

The augmented ODE state is `[n_minerals..., Q]` where `Q` is the last element.

# Fields

  - `T`: fixed temperature [K].

# Examples

```julia
cal = IsothermalCalorimeter(298.15)   # 25 °C
sol = integrate(kp, ks; calorimeter = cal)
t, Q   = cumulative_heat(sol, cal)
t, qdot = heat_flow(sol, cal)
```
"""
struct IsothermalCalorimeter{T <: Real} <: AbstractCalorimeter
    T::T
end

"""
    IsothermalCalorimeter(T_K::Real) -> IsothermalCalorimeter

Construct an [`IsothermalCalorimeter`](@ref) at temperature `T_K` [K].
Strips units if a `Quantity` is supplied.
"""
IsothermalCalorimeter(T_K) = IsothermalCalorimeter{Float64}(Float64(ustrip(T_K)))

"""
    n_extra_states(::IsothermalCalorimeter) -> Int

Number of extra ODE states added by this calorimeter (1: cumulative heat Q).
"""
n_extra_states(::IsothermalCalorimeter) = 1

"""
    extend_u0(u0::AbstractVector, ::IsothermalCalorimeter) -> Vector

Append Q₀ = 0 to the kinetic mole vector `u0`.
"""
function extend_u0(u0::AbstractVector, ::IsothermalCalorimeter)
    return vcat(u0, zero(eltype(u0)))
end

"""
    extend_ode!(du, u, p, n_kin, cal::IsothermalCalorimeter)

Append `dQ/dt = q̇` to the ODE right-hand side.

Called by `KineticsOrdinaryDiffEqExt` after the mineral ODE step.
`n_kin` is the number of kinetic species (index offset into `u` for extra states).
"""
function extend_ode!(du, ::Any, p, n_kin::Int, cal::IsothermalCalorimeter)
    # rates = -du[1:n_kin] (sign: dissolution lowers mineral → r = -du[i])
    rates = [-du[j] for j in 1:n_kin]
    qdot = heat_rate(p.kin_rxns, rates, cal.T)
    du[n_kin + 1] = qdot
    return nothing
end

# ── SemiAdiabaticCalorimeter ──────────────────────────────────────────────────

"""
    struct SemiAdiabaticCalorimeter{T<:Real, F} <: AbstractCalorimeter

Semi-adiabatic calorimeter: temperature evolves according to

```math
\\frac{dT}{dt} = \\frac{\\dot{q}(t) - \\varphi(T(t) - T_{\\rm env})}{C_p}
```

where:
- `q̇(t) = Σᵢ rᵢ × ΔHᵣ,ᵢ(T)` [W] is the instantaneous heat-generation rate,
- `φ(ΔT)` [W] is a user-supplied **heat-loss function** (e.g. linear, quadratic),
- `Cp` [J/K] is the total heat capacity of the vessel + sample,
- `T_env` [K] is the ambient temperature.

The augmented ODE state is `[n_minerals..., T]` where `T` is the last element.

The heat-loss function `φ` maps temperature difference `ΔT = T - T_env` [K] to
heat-loss rate [W]. Common choices:

| Law | Expression |
|-----|-----------|
| Linear (Newton) | `ΔT -> L * ΔT` |
| Quadratic (Lavergne et al. 2018) | `ΔT -> a*ΔT + b*ΔT^2` |
| Arbitrary | any callable `φ(ΔT)` |

# Fields

  - `T0`: initial temperature [K].
  - `T_env`: ambient temperature [K].
  - `Cp`: total heat capacity [J/K].
  - `heat_loss`: callable `ΔT::Real -> Real [W]` — heat-loss function.

# Examples

```julia
# Linear heat loss (Newton's law of cooling, L = 0.5 W/K)
cal = SemiAdiabaticCalorimeter(; T0=293.15, T_env=293.15, Cp=4000.0,
                                  heat_loss = ΔT -> 0.5 * ΔT)

# Quadratic heat loss (Lavergne et al. 2018: a [W/K], b [W/K²])
cal = SemiAdiabaticCalorimeter(; T0=293.15, T_env=293.15, Cp=4000.0,
                                  heat_loss = ΔT -> 0.48*ΔT + 0.002*ΔT^2)

# Convenience constructor with linear L keyword (backward-compatible)
cal = SemiAdiabaticCalorimeter(; T0=293.15, T_env=293.15, L=0.5, Cp=4000.0)

sol = integrate(kp, ks; calorimeter = cal)
t, T_vec = temperature_profile(sol, cal)
t, qdot  = heat_flow(sol, cal)
```

# References

  - Lavergne, F., Ben Fraj, A., Bayane, I. & Barthélémy, J.-F. (2018).
    *Estimating the mechanical properties of hydrating blended cementitious materials:
    An investigation based on micromechanics.* Cement and Concrete Research **104**, 37–60.
  - Lerch, R. (2007). Calorimetry of cement hydration. *Cem. Concr. Res.*
"""
struct SemiAdiabaticCalorimeter{T <: Real, F} <: AbstractCalorimeter
    T0::T
    T_env::T
    Cp::T
    heat_loss::F   # callable: ΔT [K] → heat loss [W]
end

"""
    SemiAdiabaticCalorimeter(; T0, T_env, Cp, heat_loss=nothing, L=nothing)

Keyword constructor. Exactly one of `heat_loss` or `L` must be provided:

  - `heat_loss`: any callable `ΔT -> [W]` (general heat-loss function).
  - `L`: linear heat-loss coefficient [W/K] (shorthand for `heat_loss = ΔT -> L*ΔT`).

# Examples

```julia
# Linear heat loss via L shorthand
cal = SemiAdiabaticCalorimeter(; T0=293.15, T_env=293.15, Cp=4000.0, L=0.5)

# Quadratic heat loss (Lavergne et al. 2018)
cal = SemiAdiabaticCalorimeter(;
    T0        = 293.15,
    T_env     = 293.15,
    Cp        = 4000.0,
    heat_loss = ΔT -> 0.48*ΔT + 0.002*ΔT^2,
)
```
"""
function SemiAdiabaticCalorimeter(;
        T0::Real,
        T_env::Real,
        Cp::Real,
        heat_loss = nothing,
        L = nothing,
    )
    hl = if !isnothing(heat_loss)
        heat_loss
    elseif !isnothing(L)
        L_f = float(L)
        ΔT -> L_f * ΔT
    else
        throw(ArgumentError("SemiAdiabaticCalorimeter requires either `heat_loss` or `L`"))
    end
    T0_f, T_env_f, Cp_f = promote(float(T0), float(T_env), float(Cp))
    return SemiAdiabaticCalorimeter{typeof(T0_f), typeof(hl)}(T0_f, T_env_f, Cp_f, hl)
end

"""
    n_extra_states(::SemiAdiabaticCalorimeter) -> Int

Number of extra ODE states added by this calorimeter (1: temperature T).
"""
n_extra_states(::SemiAdiabaticCalorimeter) = 1

"""
    extend_u0(u0::AbstractVector, cal::SemiAdiabaticCalorimeter) -> Vector

Append T₀ to the kinetic mole vector `u0`.
"""
function extend_u0(u0::AbstractVector, cal::SemiAdiabaticCalorimeter)
    return vcat(u0, eltype(u0)(cal.T0))
end

"""
    extend_ode!(du, u, p, n_kin, cal::SemiAdiabaticCalorimeter)

Append `dT/dt = (q̇ - φ(T - T_env)) / Cp` to the ODE right-hand side,
where `φ = cal.heat_loss` is the user-supplied heat-loss function.
"""
function extend_ode!(du, u, p, n_kin::Int, cal::SemiAdiabaticCalorimeter)
    T_curr = u[n_kin + 1]
    rates = [-du[j] for j in 1:n_kin]
    qdot = heat_rate(p.kin_rxns, rates, T_curr; ϵ = p.ϵ)
    ΔT = T_curr - cal.T_env
    du[n_kin + 1] = (qdot - cal.heat_loss(ΔT)) / cal.Cp
    return nothing
end

# ── Result extraction ─────────────────────────────────────────────────────────

"""
    heat_flow(sol, cal::AbstractCalorimeter) -> (t::Vector, qdot::Vector)

Extract the instantaneous heat-generation rate q̇(t) [W] from an ODE solution.

For [`IsothermalCalorimeter`](@ref): if total-enthalpy tracking data are available
in `sol.prob.p.saved_H` (filled by `KineticsOrdinaryDiffEqExt`), q̇ is derived from
H(0) - H(t), capturing both kinetic and equilibrium heat. Otherwise falls back to
finite differences on the ODE state Q (kinetic heat only).

For [`SemiAdiabaticCalorimeter`](@ref): q̇ is re-derived from the temperature ODE.

Returns `(t_vector, qdot_vector)`.
"""
function heat_flow(sol, ::IsothermalCalorimeter)
    t, Q = cumulative_heat(sol, IsothermalCalorimeter(0.0))   # reuse logic
    qdot = similar(Q)
    qdot[1] = zero(eltype(Q))
    for i in 2:lastindex(t)
        dt = t[i] - t[i - 1]
        qdot[i] = dt > 0 ? (Q[i] - Q[i - 1]) / dt : zero(eltype(Q))
    end
    return t, qdot
end

function heat_flow(sol, cal::SemiAdiabaticCalorimeter)
    t = sol.t
    n_kin = length(sol.u[1]) - n_extra_states(cal)
    # Re-derive q̇ = Cp dT/dt + φ(T - T_env) via the calorimeter energy balance
    T_vec = [u[n_kin + 1] for u in sol.u]
    qdot = similar(T_vec)
    qdot[1] = zero(eltype(T_vec))
    for i in 2:lastindex(t)
        dt = t[i] - t[i - 1]
        dTdt = dt > 0 ? (T_vec[i] - T_vec[i - 1]) / dt : zero(eltype(T_vec))
        ΔT = T_vec[i] - cal.T_env
        qdot[i] = cal.Cp * dTdt + cal.heat_loss(ΔT)
    end
    return t, qdot
end

"""
    cumulative_heat(sol, cal::IsothermalCalorimeter) -> (t::Vector, Q::Vector)

Extract the cumulative heat Q(t) = ∫₀ᵗ q̇(τ) dτ [J] from an ODE solution.

When total-enthalpy tracking data are available (i.e. `sol.prob.p.saved_H` is
non-empty — filled at each accepted ODE step by `KineticsOrdinaryDiffEqExt`),
returns:

```math
Q(t) = H(0) - H(t)
```

where H(t) = Σᵢ nᵢ(t) ΔₐH⁰ᵢ(T) is the full-system enthalpy. This captures heat
from both kinetic reactions **and** instantaneous equilibrium re-speciation.

Falls back to the ODE state Q (kinetic heat only) when enthalpy tracking is
unavailable (no species with `:ΔₐH⁰` data, or run without `OrdinaryDiffEqExt`).
"""
function cumulative_heat(sol, ::IsothermalCalorimeter)
    p = sol.prob.p
    # Prefer total-enthalpy path (includes equilibrium heat)
    if hasproperty(p, :saved_H) && !isempty(p.saved_H)
        return p.saved_t, p.saved_H[1] .- p.saved_H
    end
    # Fallback: kinetic-only heat from ODE state Q
    n_kin = length(sol.u[1]) - 1
    Q = [u[n_kin + 1] for u in sol.u]
    return sol.t, Q
end

"""
    cumulative_heat(sol, cal::SemiAdiabaticCalorimeter) -> (t::Vector, Q::Vector)

Extract the cumulative heat Q(t) = ∫₀ᵗ q̇(τ) dτ [J] from a semi-adiabatic calorimeter
ODE solution by integrating the heat-flow rate derived from the temperature ODE.

The heat-flow rate is reconstructed via the calorimeter energy balance:
`q̇(t) = Cp dT/dt + φ(T(t) - T_env)`.
"""
function cumulative_heat(sol, cal::SemiAdiabaticCalorimeter)
    t, qdot = heat_flow(sol, cal)
    Q = similar(qdot)
    Q[1] = zero(eltype(qdot))
    for i in 2:lastindex(t)
        dt = t[i] - t[i - 1]
        Q[i] = Q[i - 1] + qdot[i] * dt
    end
    return t, Q
end

"""
    temperature_profile(sol, cal::SemiAdiabaticCalorimeter) -> (t::Vector, T::Vector)

Extract the temperature profile T(t) [K] from a semi-adiabatic calorimeter ODE solution.
"""
function temperature_profile(sol, ::SemiAdiabaticCalorimeter)
    n_extra = 1
    n_kin = length(sol.u[1]) - n_extra
    T_vec = [u[n_kin + 1] for u in sol.u]
    return sol.t, T_vec
end
