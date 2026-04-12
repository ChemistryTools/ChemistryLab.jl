using DynamicQuantities
using OrderedCollections

# ── Abstract type ──────────────────────────────────────────────────────────────

"""
    abstract type AbstractRateModel end

Base type for all kinetic rate models. A concrete subtype must be callable as:

```julia
model(; T, Ω, A_surface, lna_dict=nothing, ϵ=1e-16) -> Real
```

returning the net reaction rate [mol/s] (positive = dissolution/forward,
negative = precipitation/reverse).

All concrete subtypes must be AD-compatible (ForwardDiff-safe): no `Float64`
casts in the evaluation path.
"""
abstract type AbstractRateModel end

# ── KINETICS_RATE_MODELS / KINETICS_RATE_FACTORIES ────────────────────────────

"""
    KINETICS_RATE_MODELS

Dictionary of raw kinetic rate-constant model expressions, analogous to
[`THERMO_MODELS`](@ref).

Each entry maps a model name (`:arrhenius`, …) to a `Dict` containing:
  - `:k` — symbolic `Expr` for the rate constant as a function of variables.
  - `:vars` — list of variable symbols (e.g. `[:T]`).
  - `:units` — list of `Symbol => Quantity` pairs for parameters and variables.
  - `:output_unit` — `Quantity` representing the output unit.

At package initialisation, every entry is compiled into a `ThermoFactory`
stored in [`KINETICS_RATE_FACTORIES`](@ref).

# Example

```julia
k_acid = KINETICS_RATE_FACTORIES[:arrhenius](;
    k₀    = 5.012e-1,   # mol/(m² s) at T_ref
    Ea    = 14400.0,    # J/mol
    T_ref = 298.15,     # K
)
k_acid(; T = 310.0)   # → Float64 rate constant
```
"""
const KINETICS_RATE_MODELS = Dict{Symbol, Dict}(
    :arrhenius => Dict(
        :k => :(k₀ * exp(-Ea / R_gas * (1 / T - 1 / T_ref))),
        :vars => [:T],
        :units => [
            :T => u"K",
            :T_ref => u"K",
            :k₀ => u"mol/(m^2*s)",
            :Ea => u"J/mol",
            :R_gas => u"J/(mol*K)",
        ],
        :output_unit => u"mol/(m^2*s)",
    ),
)

"""
    KINETICS_RATE_FACTORIES

Compiled `ThermoFactory` objects for each kinetic rate model.
Populated by `__init__()` from [`KINETICS_RATE_MODELS`](@ref).

Keys are model name symbols (e.g. `:arrhenius`).
Values are `ThermoFactory` callables that return `SymbolicFunc{1}` instances.

# Usage

```julia
factory = KINETICS_RATE_FACTORIES[:arrhenius]
k = factory(; k₀=1e-5, Ea=50000.0, T_ref=298.15, R_gas=8.31446)
k(; T = 298.15)   # → 1e-5  (rate constant at reference temperature)
```
"""
const KINETICS_RATE_FACTORIES = Dict{Symbol, ThermoFactory}()

"""
    add_kinetics_rate_model(name::Symbol, dict_model::Dict)

Register a new kinetic rate-constant model in [`KINETICS_RATE_MODELS`](@ref)
and compile it into [`KINETICS_RATE_FACTORIES`](@ref).

`dict_model` must contain at minimum `:k` (expression), `:vars` (variable list),
`:units` (parameter units), and `:output_unit`.

# Example

```julia
add_kinetics_rate_model(:power_law, Dict(
    :k    => :(k₀ * (T / T_ref)^n),
    :vars => [:T],
    :units => [:T => u"K", :T_ref => u"K", :k₀ => u"mol/(m^2*s)", :n => u"1"],
    :output_unit => u"mol/(m^2*s)",
))
```
"""
function add_kinetics_rate_model(name::Symbol, dict_model::Dict)
    KINETICS_RATE_MODELS[name] = dict_model
    KINETICS_RATE_FACTORIES[name] = _build_kinetics_rate_factory(dict_model)
    return nothing
end

# Internal: compile one KINETICS_RATE_MODELS entry → ThermoFactory
function _build_kinetics_rate_factory(d::Dict)
    return ThermoFactory(
        d[:k],
        get(d, :vars, [:T]);
        units = get(d, :units, nothing),
        output_unit = get(d, :output_unit, u"1"),
    )
end

# ── arrhenius_rate_constant ────────────────────────────────────────────────────

"""
    arrhenius_rate_constant(k₀, Ea; T_ref=298.15, R_gas=8.31446261815324) -> NumericFunc

Build a temperature-dependent Arrhenius rate constant as a [`NumericFunc`](@ref):

```
k(T) = k₀ × exp(-Eₐ / R × (1/T - 1/T_ref))
```

The returned object is callable as `k(; T=...)` and fully AD-compatible
(ForwardDiff-safe: the closure captures `k₀`, `Ea`, `T_ref`, `R_gas` directly,
so dual numbers propagate correctly through all parameters).

Arithmetic between `SymbolicFunc`/`NumericFunc` objects is supported, so rate
constants can be composed with activity or surface-area functions.

# Arguments

  - `k₀`: pre-exponential factor [mol/(m² s)] at `T_ref`.
  - `Ea`: activation energy [J/mol].
  - `T_ref`: reference temperature [K] (default `298.15`).
  - `R_gas`: gas constant [J/(mol K)] (default `8.31446261815324`).

# Returns

A `NumericFunc` with variable `T` (in K) and `refs = (T = T_ref * u"K",)`.

# Examples

```jldoctest
julia> k = arrhenius_rate_constant(5.0e-4, 62000.0);

julia> isapprox(k(; T = 298.15), 5.0e-4; rtol = 1e-10)
true

julia> k(; T = 350.0) > k(; T = 298.15)   # higher T → higher k
true
```

AD-compatible through all parameters:
```julia
ForwardDiff.derivative(T  -> arrhenius_rate_constant(5e-4, 62000.0)(; T = T),  298.15)
ForwardDiff.derivative(Ea -> arrhenius_rate_constant(5e-4, Ea)(; T = 350.0),   62000.0)
ForwardDiff.derivative(k₀ -> arrhenius_rate_constant(k₀,   62000.0)(; T = 298.15), 5e-4)
```
"""
function arrhenius_rate_constant(
        k₀::Real,
        Ea::Real;
        T_ref::Real = 298.15,
        R_gas::Real = 8.31446261815324,
    )
    # NumericFunc closure: captures k₀, Ea, T_ref, R_gas → ForwardDiff-safe
    # through all parameters (ThermoFactory would bake them as constants).
    f = (T,) -> k₀ * exp(-Ea / R_gas * (1 / T - 1 / T_ref))
    refs = (T = T_ref * u"K",)
    return NumericFunc(f, (:T,), refs, u"mol/(m^2*s)")
end

# ── Saturation ratio ───────────────────────────────────────────────────────────

"""
    saturation_ratio(stoich::AbstractVector, lna::AbstractVector,
                     ΔₐG⁰overT::AbstractVector, T_K; ϵ=1e-16) -> Real

Compute the saturation ratio Ω = IAP / K for a kinetic reaction.

```
ln Ω = Σᵢ νᵢ ln aᵢ + ln K
     = Σᵢ νᵢ ln aᵢ + ΔᵣG⁰/(RT)   (note: ΔᵣG⁰/RT = Σᵢ νᵢ ΔₐG⁰ᵢ/RT for reactants→products)
```

where `stoich[i]` is the stoichiometric coefficient (positive for products,
negative for reactants), `lna[i]` is the log-activity of species `i`,
and `ΔₐG⁰overT[i]` is the dimensionless standard Gibbs energy of formation
`ΔₐG⁰ᵢ / RT` for species `i`.

# Arguments

  - `stoich`: stoichiometric coefficient vector for this reaction (length = number of species).
  - `lna`: log-activity vector (same indexing as species in system).
  - `ΔₐG⁰overT`: dimensionless standard Gibbs energies `ΔₐG⁰ᵢ/RT`.
  - `ϵ`: floor to avoid `exp` overflow when Ω → ∞.

# Returns

`Ω = exp(ln_IAP - ln_K)` where `ln_K = -ΔᵣG⁰/RT`.

AD-compatible (ForwardDiff-safe).
"""
function saturation_ratio(
        stoich::AbstractVector,
        lna::AbstractVector,
        ΔₐG⁰overT::AbstractVector;
        ϵ::Real = 1.0e-16,
    )
    # ln IAP = Σᵢ νᵢ ln aᵢ
    ln_iap = sum(stoich[i] * lna[i] for i in eachindex(stoich))
    # ln K = -ΔᵣG⁰/RT = -Σᵢ νᵢ ΔₐG⁰ᵢ/RT
    ln_K = -sum(stoich[i] * ΔₐG⁰overT[i] for i in eachindex(stoich))
    return exp(ln_iap - ln_K)
end

# ── RateModelCatalyst ─────────────────────────────────────────────────────────

"""
    struct RateModelCatalyst{T<:Real}

Describes the contribution of a catalyst species to a reaction mechanism rate.

The catalyst multiplies the base rate by `exp(n * ln aᵢ) = aᵢ^n`, where
`aᵢ` is the activity of the catalyst species.

# Fields

  - `species`: PHREEQC-format formula string of the catalyst species (e.g. `"H+"`, `"OH-"`).
  - `n`: power exponent (dimensionless).

# Examples

```julia
acid_catalyst   = RateModelCatalyst("H+",  0.5)    # ∝ a(H+)^0.5
base_catalyst   = RateModelCatalyst("OH-", 0.5)    # ∝ a(OH-)^0.5
co2_catalyst    = RateModelCatalyst("CO2", 1.0)    # ∝ a(CO2)
```
"""
struct RateModelCatalyst{T <: Real}
    species::String
    n::T
end

# ── RateMechanism ─────────────────────────────────────────────────────────────

"""
    struct RateMechanism{F<:AbstractFunc, T<:Real}

A single kinetic mechanism (acid/neutral/base/…) contributing to the overall
mineral dissolution or precipitation rate.

The mechanism rate is:
```
r_mech = k(T) × [Π_catalysts aᵢ^nᵢ] × sign(1 - Ω) × |1 - Ω^p|^q
```

# Fields

  - `k`: rate constant as `AbstractFunc` (typically `SymbolicFunc{1}` from
    [`arrhenius_rate_constant`](@ref)). Called as `k(; T=...)`.
  - `p`: saturation exponent `p` in `(1 - Ω^p)^q`. Default 1.0.
  - `q`: outer exponent `q`. Default 1.0.
  - `catalysts`: vector of [`RateModelCatalyst`](@ref) (may be empty).

# Examples

```julia
k_acid = arrhenius_rate_constant(5.012e-1, 14400.0)
mech   = RateMechanism(k_acid, 1.0, 1.0, [RateModelCatalyst("H+", 1.0)])
```
"""
struct RateMechanism{F <: AbstractFunc, T <: Real}
    k::F
    p::T
    q::T
    catalysts::Vector{RateModelCatalyst{T}}
end

"""
    RateMechanism(k::AbstractFunc, p::Real, q::Real) -> RateMechanism

Construct a [`RateMechanism`](@ref) with no catalyst contributions.
"""
function RateMechanism(k::AbstractFunc, p::Real, q::Real)
    T = typeof(promote(p, q)[1])
    return RateMechanism{typeof(k), T}(k, T(p), T(q), RateModelCatalyst{T}[])
end

# ── TransitionStateRateModel ──────────────────────────────────────────────────

"""
    struct TransitionStateRateModel{T<:Real} <: AbstractRateModel

Multi-mechanism mineral dissolution / precipitation rate model based on
transition-state theory (Palandri & Kharaka 2004).

The net rate is:
```
r = A_surface × Σ_mechanisms r_mech
  = A_surface × Σ_m [ k_m(T) × Π_cat(aᵢ^nᵢ) × sign(1-Ω) × |1-Ω^p|^q ]
```

where Ω = IAP/K is the saturation ratio.  Ω < 1 → dissolution (r > 0);
Ω > 1 → precipitation (r < 0).

# Fields

  - `mechanisms`: vector of [`RateMechanism`](@ref) (acid, neutral, base, …).

# References

  - Palandri, J.L. & Kharaka, Y.K. (2004). USGS Open-File Report 2004-1068.
  - Leal, A.M.M. et al. (2017). Pure Appl. Chem. 89, 597–643.
    https://doi.org/10.1515/pac-2016-1107

# Examples

```julia
k_acid    = arrhenius_rate_constant(5.012e-1, 14400.0)
k_neutral = arrhenius_rate_constant(1.549e-6, 23500.0)

model = TransitionStateRateModel([
    RateMechanism(k_acid,    1.0, 1.0, [RateModelCatalyst("H+", 1.0)]),
    RateMechanism(k_neutral, 1.0, 1.0),
])
```
"""
struct TransitionStateRateModel{T <: Real} <: AbstractRateModel
    mechanisms::Vector{<:RateMechanism}
end

"""
    TransitionStateRateModel(mechanisms::AbstractVector) -> TransitionStateRateModel

Construct a [`TransitionStateRateModel`](@ref) from a vector of mechanisms.
"""
function TransitionStateRateModel(mechanisms::AbstractVector{<:RateMechanism})
    T = isempty(mechanisms) ? Float64 : typeof(mechanisms[1].p)
    return TransitionStateRateModel{T}(collect(mechanisms))
end

"""
    (model::TransitionStateRateModel)(; T, Ω, A_surface, lna_dict=nothing, ϵ=1e-16) -> Real

Evaluate the net reaction rate [mol/s].

# Arguments

  - `T`: temperature in K (plain number; strips units internally).
  - `Ω`: saturation ratio IAP/K (dimensionless; see [`saturation_ratio`](@ref)).
  - `A_surface`: reactive surface area [m²].
  - `lna_dict`: `Dict{String, Real}` mapping PHREEQC species formulas to their
    log-activities. Required only when any mechanism has catalysts. Optional.
  - `ϵ`: regularization floor (default `1e-16`) to smooth the rate near Ω = 1.

# Returns

Net rate in mol/s (positive = dissolution, negative = precipitation).
AD-compatible: no `Float64` casts in the evaluation path.
"""
function (model::TransitionStateRateModel)(;
        T,
        Ω,
        A_surface,
        lna_dict::Union{Nothing, AbstractDict} = nothing,
        ϵ::Real = 1.0e-16,
        kwargs...,   # accept n_current, n_initial, etc. from generic ODE dispatch
    )
    T_val = ustrip(T)    # strip units if Quantity; plain numbers pass through
    Ω_val = ustrip(Ω)

    r = zero(promote_type(typeof(T_val), typeof(Ω_val), typeof(A_surface)))

    for mech in model.mechanisms
        k_val = mech.k(; T = T_val)

        # Catalyst product: Π aᵢ^nᵢ = exp(Σ nᵢ ln aᵢ)
        cat_term = one(typeof(r))
        if !isnothing(lna_dict) && !isempty(mech.catalysts)
            for cat in mech.catalysts
                ln_a = get(lna_dict, cat.species, zero(typeof(r)))
                cat_term *= exp(cat.n * ln_a)
            end
        end

        # Saturation driving force: sign(1-Ω) * |1 - Ω^p|^q
        # Smooth form: avoids branching on Dual values
        Ωp = Ω_val^mech.p
        diff = one(typeof(r)) - Ωp
        # sign-preserving smooth formulation for AD:
        sat_term = diff * (diff^2 + ϵ)^((mech.q - one(typeof(r))) / 2)

        r = r + k_val * cat_term * sat_term
    end

    return A_surface * r
end

# ── FirstOrderRateModel ───────────────────────────────────────────────────────

"""
    struct FirstOrderRateModel{F<:AbstractFunc, T<:Real} <: AbstractRateModel

Simple first-order dissolution/precipitation rate model:

```
r = A_surface × k(T) × sign(1 - Ω) × |1 - Ω^p|^q
```

Useful as a minimal test case or for empirical fits.

# Fields

  - `k`: rate constant as `AbstractFunc` (e.g. from [`arrhenius_rate_constant`](@ref)).
  - `p`: saturation exponent `p` in `|1 - Ω^p|^q`. Default 1.0.
  - `q`: outer exponent `q`. Default 1.0.

# Examples

```julia
k = arrhenius_rate_constant(1e-7, 40000.0)
rm = FirstOrderRateModel(k)
rm(; T = 298.15, Ω = 0.5, A_surface = 1.0)   # dissolution rate [mol/s]
```
"""
struct FirstOrderRateModel{F <: AbstractFunc, T <: Real} <: AbstractRateModel
    k::F
    p::T
    q::T
end

"""
    FirstOrderRateModel(k::AbstractFunc; p=1.0, q=1.0) -> FirstOrderRateModel

Construct a [`FirstOrderRateModel`](@ref) with optional saturation exponents.
"""
FirstOrderRateModel(k::AbstractFunc; p::Real = 1.0, q::Real = 1.0) =
    FirstOrderRateModel{typeof(k), typeof(float(p))}(k, float(p), float(q))

"""
    (model::FirstOrderRateModel)(; T, Ω, A_surface, ϵ=1e-16, kwargs...) -> Real

Evaluate the first-order rate. AD-compatible.
Extra keyword arguments (e.g. `lna_dict`) are accepted and silently ignored
for interface compatibility with [`TransitionStateRateModel`](@ref).
"""
function (model::FirstOrderRateModel)(;
        T,
        Ω,
        A_surface,
        ϵ::Real = 1.0e-16,
        kwargs...,
    )
    T_val = ustrip(T)
    Ω_val = ustrip(Ω)
    k_val = model.k(; T = T_val)
    Ωp = Ω_val^model.p
    diff = one(typeof(k_val)) - Ωp
    sat = diff * (diff^2 + ϵ)^((model.q - one(typeof(k_val))) / 2)
    return A_surface * k_val * sat
end

# ── ParrotKillohRateModel ────────────────────────────────────────────────────

"""
    struct ParrotKillohRateModel{T<:Real} <: AbstractRateModel

Parrot & Killoh (1984) kinetic model for cement clinker hydration.

Three competing mechanisms determine the hydration rate:

```math
\\frac{d\\alpha}{dt} = A_T \\cdot \\min(\\max(r_{\\rm NG},\\, r_{\\rm I}),\\, r_{\\rm D})
```

where ``\\alpha`` is the degree of hydration (0 to ``\\alpha_{\\rm max}``),
``\\xi = \\alpha / \\alpha_{\\rm max}`` is the normalized degree of hydration,
and ``A_T = \\exp(-E_a/R \\cdot (1/T - 1/T_{\\rm ref}))`` is the Arrhenius factor.
The three mechanisms are:

| Mechanism | Formula |
|-----------|---------|
| Nucleation–growth | ``r_{\\rm NG} = (K_1/N_1)(1-\\xi)^{N_1} / (1 + B\\xi^{N_3})`` |
| Interaction | ``r_{\\rm I} = K_2(1-\\xi)^{N_2}`` |
| Diffusion | ``r_{\\rm D} = 3K_3(1-\\xi)^{2/3} / (N_3(1-(1-\\xi)^{1/3}))`` |

Rate constants ``K_1, K_2, K_3`` are in d⁻¹ (per day); conversion to s⁻¹ is
applied internally.

The model requires the current and initial mineral moles (keywords `n_current`,
`n_initial`) and ignores `Ω`, `A_surface`, and `lna_dict` (it is not
surface-area controlled and does not use saturation ratio).

# Fields

  - `K₁`, `N₁`: nucleation–growth rate [d⁻¹] and exponent.
  - `K₂`, `N₂`: interaction rate [d⁻¹] and exponent.
  - `K₃`, `N₃`: diffusion rate [d⁻¹] and geometry factor.
  - `B`: geometry exponent for the nucleation–growth denominator.
  - `Ea`: activation energy [J/mol].
  - `T_ref`: reference temperature [K] for the Arrhenius correction (default 293.15 K).
  - `α_max`: maximum reachable degree of hydration (0 < α_max ≤ 1).

# Examples

```julia
using ChemistryLab

# Use predefined C3S parameters (Parrot & Killoh 1984)
rm = PK_PARAMS_C3S

# Evaluate at 10 % hydration, 25 °C, for 1 mol initial C3S
rm(; T = 298.15, n_current = 0.9, n_initial = 1.0)   # mol/s
```

See also: [`KineticReaction`](@ref), [`PK_PARAMS_C3S`](@ref),
[`PK_PARAMS_C2S`](@ref), [`PK_PARAMS_C3A`](@ref), [`PK_PARAMS_C4AF`](@ref).
"""
struct ParrotKillohRateModel{T <: Real} <: AbstractRateModel
    K₁::T    # nucleation–growth rate constant [d⁻¹]
    N₁::T    # nucleation–growth exponent
    K₂::T    # interaction rate constant [d⁻¹]
    N₂::T    # interaction exponent
    K₃::T    # diffusion rate constant [d⁻¹]
    N₃::T    # diffusion geometry factor
    B::T     # geometry exponent for NG denominator
    Ea::T    # activation energy [J/mol]
    T_ref::T # reference temperature [K]
    α_max::T # maximum degree of hydration
end

"""
    ParrotKillohRateModel(K₁, N₁, K₂, N₂, K₃, N₃, B, Ea;
                           T_ref=293.15, α_max=1.0) -> ParrotKillohRateModel

Construct a [`ParrotKillohRateModel`](@ref).

# Arguments

  - `K₁`, `N₁`: nucleation–growth rate [d⁻¹] and exponent.
  - `K₂`, `N₂`: interaction rate [d⁻¹] and exponent.
  - `K₃`, `N₃`: diffusion rate [d⁻¹] and geometry factor.
  - `B`: geometry exponent for the NG denominator.
  - `Ea`: activation energy [J/mol].
  - `T_ref`: reference temperature [K] (default: 293.15 = 20 °C).
  - `α_max`: maximum degree of hydration (default: 1.0).
"""
function ParrotKillohRateModel(
        K₁::Real, N₁::Real, K₂::Real, N₂::Real,
        K₃::Real, N₃::Real, B::Real, Ea::Real;
        T_ref::Real = 293.15,
        α_max::Real = 1.0,
    )
    T = promote_type(
        typeof(float(K₁)), typeof(float(N₁)),
        typeof(float(K₂)), typeof(float(N₂)),
        typeof(float(K₃)), typeof(float(N₃)),
        typeof(float(B)), typeof(float(Ea)),
        typeof(float(T_ref)), typeof(float(α_max)),
    )
    return ParrotKillohRateModel{T}(
        T(K₁), T(N₁), T(K₂), T(N₂),
        T(K₃), T(N₃), T(B), T(Ea),
        T(T_ref), T(α_max),
    )
end

"""
    (model::ParrotKillohRateModel)(; T, n_current, n_initial, ϵ=1e-10, kwargs...) -> Real

Evaluate the Parrot–Killoh hydration rate [mol/s].

# Arguments

  - `T`: temperature [K].
  - `n_current`: current moles of the unreacted clinker mineral.
  - `n_initial`: initial moles of the same mineral at t = 0.
  - `ϵ`: regularisation floor for the diffusion denominator (default `1e-10`).

Extra keyword arguments (`Ω`, `A_surface`, `lna_dict`, …) are accepted and
silently ignored for interface compatibility with other rate models.

# Returns

Rate [mol/s] of clinker consumption (positive = dissolution).
AD-compatible: no `Float64` casts in the evaluation path.
"""
function (model::ParrotKillohRateModel)(;
        T,
        n_current,
        n_initial,
        ϵ::Real = 1.0e-10,
        kwargs...,   # accept Ω, A_surface, lna_dict, etc. — ignored
    )
    T_val = ustrip(T)
    n0 = ustrip(n_initial)
    nc = ustrip(n_current)

    # degree of hydration α ∈ [0, α_max)
    α = min(max(one(T_val) - nc / n0, zero(T_val)), model.α_max - ϵ)
    ξ = α / model.α_max

    # Arrhenius temperature correction
    R_gas = 8.314462   # J/(mol·K)
    Aₜ = exp(-model.Ea / R_gas * (one(T_val) / T_val - one(model.T_ref) / model.T_ref))

    one_m_ξ = one(ξ) - ξ

    # r_NG: nucleation–growth  [s⁻¹]
    r_NG = (model.K₁ / model.N₁) * one_m_ξ^model.N₁ / (one(ξ) + model.B * ξ^model.N₃) / 86400

    # r_I: interaction (dissolution–precipitation)  [s⁻¹]
    r_I = model.K₂ * one_m_ξ^model.N₂ / 86400

    # r_D: diffusion  [s⁻¹]  (denominator clamped to avoid 0/0 at α = 0)
    denom_D = max(one(ξ) - one_m_ξ^(one(ξ) / 3), ϵ)
    r_D = 3 * model.K₃ * one_m_ξ^(2 * one(ξ) / 3) / (model.N₃ * denom_D) / 86400

    # dα/dt [s⁻¹] at temperature T
    dα_dt = Aₜ * min(max(r_NG, r_I), r_D)

    # Convert to mol/s: dn/dt = -n₀ · dα/dt  →  dissolution rate = n₀ · dα/dt
    return n0 * dα_dt
end

# ── Predefined Parrot & Killoh (1984) parameters ────────────────────────────

"""
    PK_PARAMS_C3S :: ParrotKillohRateModel

Parrot & Killoh (1984) parameters for alite (C₃S = Ca₃SiO₅).

Activation energy from Schindler & Folliard (2005).
Reference temperature: 293.15 K (20 °C).
"""
const PK_PARAMS_C3S = ParrotKillohRateModel(
    1.5, 3.3, 0.018, 2.5, 0.0024, 4.0, 0.5, 41_570.0;
    T_ref = 293.15, α_max = 1.0,
)

"""
    PK_PARAMS_C2S :: ParrotKillohRateModel

Parrot & Killoh (1984) parameters for belite (C₂S = Ca₂SiO₄).

Activation energy from Schindler & Folliard (2005).
Reference temperature: 293.15 K (20 °C).
"""
const PK_PARAMS_C2S = ParrotKillohRateModel(
    0.95, 0.5, 0.0005, 2.5, 0.0024, 4.0, 0.2, 43_670.0;
    T_ref = 293.15, α_max = 1.0,
)

"""
    PK_PARAMS_C3A :: ParrotKillohRateModel

Parrot & Killoh (1984) parameters for tricalcium aluminate (C₃A = Ca₃Al₂O₆)
in the presence of sulfate (gypsum), corresponding to ettringite formation.

Activation energy from Schindler & Folliard (2005).
Reference temperature: 293.15 K (20 °C).
"""
const PK_PARAMS_C3A = ParrotKillohRateModel(
    0.082, 0.87, 0.00024, 2.0, 0.0024, 4.0, 0.04, 54_040.0;
    T_ref = 293.15, α_max = 1.0,
)

"""
    PK_PARAMS_C4AF :: ParrotKillohRateModel

Parrot & Killoh (1984) parameters for tetracalcium aluminoferrite
(C₄AF = Ca₄Al₂Fe₂O₁₀).

Activation energy from Schindler & Folliard (2005).
Reference temperature: 293.15 K (20 °C).
"""
const PK_PARAMS_C4AF = ParrotKillohRateModel(
    0.165, 3.7, 0.0015, 2.5, 0.0024, 4.0, 0.5, 34_420.0;
    T_ref = 293.15, α_max = 1.0,
)
