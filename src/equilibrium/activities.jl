# ── Abstract activity model ───────────────────────────────────────────────────

"""
    abstract type AbstractActivityModel end

Base type for all activity models. Each concrete subtype must implement
`activity_model(cs::ChemicalSystem, model::AbstractActivityModel)`
which returns a closure `a(n, p) -> Vector{Float64}` of log-activities.
"""
abstract type AbstractActivityModel end

# ── Concrete models ───────────────────────────────────────────────────────────

"""
    struct DiluteSolutionModel <: AbstractActivityModel

Ideal dilute solution model:
- Solvent:  Raoult's law  — `ln a = ln(x_solvent)`
- Solutes:  Henry's law   — `ln a = ln(c_i / c°)` where `c° = 1 mol/L`
- Crystals: pure solid    — `ln a = 0`
- Gas:      ideal mixture — `ln a = ln(x_i)`
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

All quantities are dimensionless — units are stripped at construction time.
"""
function activity_model(cs::ChemicalSystem, ::DiluteSolutionModel)

    idx_solvent  = only(cs.idx_solvent)
    idx_solutes  = cs.idx_solutes
    idx_crystal  = cs.idx_crystal
    idx_gas      = cs.idx_gas

    M_solvent    = ustrip(us"kg/mol", cs.species[idx_solvent][:M])
    ρ_solvent    = 1.0
    c_solvent    = ρ_solvent / M_solvent
    ln_c_solvent = log(c_solvent)

    function lna(n::AbstractVector, p)
        ϵ  = p.ϵ
        _n = max.(n, ϵ)     # ϵ::Float64 — promotion vers Dual automatique si n est Dual

        out = zeros(eltype(_n), length(_n))

        n_aqueous = _n[idx_solvent] + sum((_n[i] for i in idx_solutes); init=0.0)
        if !iszero(ustrip(n_aqueous))
            out[idx_solvent] = log(_n[idx_solvent] / n_aqueous)
            for i in idx_solutes
                out[i] = log(_n[i] / _n[idx_solvent]) + ln_c_solvent
            end
        end

        n_gas = sum((_n[i] for i in idx_gas); init=0.0)
        if !iszero(ustrip(n_gas))
            for i in idx_gas
                out[i] = log(_n[i] / n_gas)
            end
        end

        return out
    end

    return lna
end

# ── HKF (B-dot / Extended Debye-Hückel) activity model ───────────────────────

# ── Ion-size database ─────────────────────────────────────────────────────────

"""
    _rej_to_phreeqc(formula) -> String

Convert a species formula from `rej.dat` notation (e.g. `"Ca++"`, `"SO4--"`)
to ChemistryLab PHREEQC notation (e.g. `"Ca+2"`, `"SO4-2"`).
"""
function _rej_to_phreeqc(formula::AbstractString)
    s = replace(formula, "+++" => "+3")
    s = replace(s, "++" => "+2")
    s = replace(s, "--" => "-2")
    return s
end

"""
    _parse_rej_dat(filename) -> Dict{String, Float64}

Parse a `rej.dat` file and return a `Dict` mapping PHREEQC species formulas to
effective electrostatic radii åⱼ [Å].  Lines not matching `<formula>  <value>`
are silently skipped.
"""
function _parse_rej_dat(filename::AbstractString)
    dict = Dict{String, Float64}()
    isfile(filename) || return dict
    for line in eachline(filename)
        m = match(r"^(\S+)\s+([\d.]+)\s*$", line)
        m === nothing && continue
        dict[_rej_to_phreeqc(m[1])] = parse(Float64, m[2])
    end
    return dict
end

"""
    REJ_HKF :: Dict{String, Float64}

Effective electrostatic radii åⱼ [Å] from Table 3 in Helgeson, Kirkham &
Flowers (1981), loaded from `data/rej.dat` at package load time.  Keys are
PHREEQC-format species formulas (e.g. `"Ca+2"`, `"SO4-2"`).

Used automatically by [`activity_model`](@ref) for [`HKFActivityModel`](@ref).
"""
const REJ_HKF = _parse_rej_dat(joinpath(@__DIR__, "..", "..", "data", "rej.dat"))

"""
    REJ_CHARGE_DEFAULT :: Dict{Int, Float64}

Charge-based fallback ion-size parameters åⱼ [Å] for species absent from
[`REJ_HKF`](@ref).  Values from Table 1 of the ToughReact V2 manual
(Xu et al. 2011), itself derived from Helgeson et al. (1981).

| Charge | å (Å) | Source                                        |
|:------:|:-----:|:----------------------------------------------|
| −1     | 1.81  | Cl⁻ value                                     |
| −2     | 3.00  | Mean of CO₃²⁻ and SO₄²⁻                       |
| −3     | 4.20  | Straight-line fit with charge                 |
| +1     | 2.31  | NH₄⁺ value                                    |
| +2     | 2.80  | Mean of +2 species in HKF Table 3             |
| +3     | 3.60  | Mean of +3 species in HKF Table 3             |
| +4     | 4.50  | HKF Eq. 142 + CRC mean crystallographic radii |
"""
const REJ_CHARGE_DEFAULT = Dict{Int, Float64}(
    -1 => 1.81,
    -2 => 3.00,
    -3 => 4.20,
     1 => 2.31,
     2 => 2.80,
     3 => 3.60,
     4 => 4.50,
)

# ─────────────────────────────────────────────────────────────────────────────

"""
    _hkf_sigma(x) -> Real

Auxiliary function for the Debye-Hückel osmotic-coefficient integral:

```math
\\sigma(x) = \\frac{3}{x^3}\\!\\left[x - 2\\ln(1+x) - \\frac{1}{1+x} + 1\\right]
```

with ``\\sigma(0) = 1``. For ``|x| < 10^{-3}`` the Taylor expansion
``1 - \\tfrac{3}{2}x + \\tfrac{9}{5}x^2`` is used to avoid catastrophic cancellation.

# Reference

Helgeson, H.C. (1969), *Am. J. Sci.* **267**, 729–804.
"""
@inline function _hkf_sigma(x::T) where {T<:Real}
    abs(x) < 1e-3 && return one(T) - T(3 / 2) * x + T(9 / 5) * x^2
    return (3 / x^3) * (x - 2 * log(1 + x) - 1 / (1 + x) + 1)
end

"""
    struct HKFActivityModel <: AbstractActivityModel

Extended Debye-Hückel (B-dot / Helgeson-Kirkham-Flowers) activity model for
aqueous solutions.

## Activity coefficients

For **charged** aqueous solutes (``z_i \\neq 0``):

```math
\\log_{10}\\gamma_i = -\\frac{A\\,z_i^2\\sqrt{I}}{1 + B\\,\\dot{a}_i\\sqrt{I}} + \\dot{B}\\,I
```

For **neutral** aqueous solutes (``z_i = 0``):

```math
\\log_{10}\\gamma_i = K_n\\,I
```

Log-activity: ``\\ln a_i = \\ln 10 \\cdot \\log_{10}(\\gamma_i\\,m_i)``
(molality standard state ``m^\\circ = 1`` mol/kg).

## Water activity

Computed via the osmotic coefficient ``\\varphi``:

```math
\\ln a_{\\mathrm{w}} = -M_{\\mathrm{w}}\\,\\sum_j m_j\\,\\varphi
```

```math
\\varphi = 1
  - \\frac{A\\ln 10}{3}\\cdot\\frac{2I}{\\sum_j m_j}\\cdot\\sqrt{I}\\cdot
    \\sigma\\!\\left(B\\,\\bar{a}\\sqrt{I}\\right)
  + \\frac{\\dot{B}\\ln 10}{2}\\,I
```

where ``\\sigma(x)`` is defined in `_hkf_sigma` and
``\\bar{a} = (\\sum_j m_j z_j^2\\,\\dot{a}_j) / (\\sum_j m_j z_j^2)`` is the
charge-squared-weighted mean ion-size parameter.

## Parameters

| Field        | Description                                      | Default (25 °C, 1 bar) |
|:-------------|:-------------------------------------------------|:----------------------:|
| `A`          | Debye-Hückel A [(kg/mol)^(1/2)]                  | 0.5114                 |
| `B`          | Debye-Hückel B [Å⁻¹ (kg/mol)^(1/2)]             | 0.3288                 |
| `Bdot`       | Extended term Ḃ [kg/mol]                         | 0.041                  |
| `Kn`         | Neutral-species salting coefficient [kg/mol]     | 0.1                    |
| `å_default`  | Fallback ion-size parameter [Å]                  | 3.72                   |

Ion-size parameters åᵢ are resolved in the following priority order:
1. `sp[:å]` [Å] if explicitly set in species properties.
2. [`REJ_HKF`](@ref) — Table 3 of Helgeson et al. (1981), keyed by PHREEQC formula.
3. [`REJ_CHARGE_DEFAULT`](@ref) — charge-based fallback (ToughReact table).
4. `å_default` as a last resort.

## Validity

Best suited for ionic strengths ``I \\lesssim 1`` mol/kg. For higher
concentrations consider the Pitzer model.

## References

- Helgeson, H.C. (1969). *Am. J. Sci.* **267**, 729–804.
- Helgeson, H.C., Kirkham, D.H. & Flowers, G.C. (1981). *Am. J. Sci.* **281**, 1249–1516.
- Parkhurst, D.L. & Appelo, C.A.J. (2013). *USGS Techniques Methods*, Book 6, A43.
"""
struct HKFActivityModel <: AbstractActivityModel
    A::Float64          # Debye-Hückel A  [(kg/mol)^(1/2)]
    B::Float64          # Debye-Hückel B  [Å⁻¹ (kg/mol)^(1/2)]
    Bdot::Float64       # B-dot           [kg/mol]
    Kn::Float64         # neutral salting [kg/mol]
    å_default::Float64  # fallback å      [Å]
end

"""
    HKFActivityModel(; A=0.5114, B=0.3288, Bdot=0.041, Kn=0.1, å_default=3.72)

Keyword constructor with default parameters at 25 °C, 1 bar (H₂O solvent).
"""
function HKFActivityModel(;
    A         = 0.5114,
    B         = 0.3288,
    Bdot      = 0.041,
    Kn        = 0.1,
    å_default = 3.72,
)
    return HKFActivityModel(A, B, Bdot, Kn, å_default)
end

"""
    activity_model(cs::ChemicalSystem, model::HKFActivityModel) -> Function

Return a closure `lna(n, p) -> Vector` of log-activities for the HKF
extended Debye-Hückel (B-dot) model.

Conventions:
- Aqueous solutes : `ln aᵢ = ln(γᵢ mᵢ)` (molality standard state, m° = 1 mol/kg)
- Solvent (water) : `ln a_w = −Mw Σmⱼ φ`  (osmotic coefficient)
- Crystals        : `ln a = 0`  (pure-solid reference)
- Gas mixture     : `ln aᵢ = ln(xᵢ)`  (ideal-gas reference)

Charges are read from `sp.formula.charge`. Ion-size parameters åᵢ are resolved
at construction time (priority: `sp[:å]` → [`REJ_HKF`](@ref) →
[`REJ_CHARGE_DEFAULT`](@ref) → `model.å_default`).
"""
function activity_model(cs::ChemicalSystem, model::HKFActivityModel)

    idx_solvent = only(cs.idx_solvent)
    idx_solutes = cs.idx_solutes
    idx_gas     = cs.idx_gas

    A    = model.A
    B    = model.B
    Bdot = model.Bdot
    Kn   = model.Kn
    ln10 = log(10)

    M_w = ustrip(us"kg/mol", cs.species[idx_solvent][:M])  # Mw [kg/mol]

    # Precompute per-solute charges and å values (captured once at construction)
    # Priority: sp[:å] > REJ_HKF (rej.dat) > REJ_CHARGE_DEFAULT (by charge) > å_default
    charges = Int[cs.species[i].formula.charge for i in idx_solutes]
    å_vec = Float64[
        let sp = cs.species[i], z = Int(sp.formula.charge)
            if haskey(sp, :å)
                Float64(sp[:å])
            elseif haskey(REJ_HKF, sp.formula.phreeqc)
                REJ_HKF[sp.formula.phreeqc]
            elseif haskey(REJ_CHARGE_DEFAULT, z)
                REJ_CHARGE_DEFAULT[z]
            else
                model.å_default
            end
        end
        for i in idx_solutes
    ]

    function lna(n::AbstractVector, p)
        ϵ  = p.ϵ
        _n = max.(n, ϵ)
        out = zeros(eltype(_n), length(_n))

        # ── Molalities [mol/kg] ────────────────────────────────────────────────
        kgw = _n[idx_solvent] * M_w                                    # kg of water
        m   = [_n[idx_solutes[k]] / kgw for k in eachindex(idx_solutes)]

        # ── Ionic strength [mol/kg] ────────────────────────────────────────────
        I     = 0.5 * sum(m[k] * charges[k]^2 for k in eachindex(idx_solutes); init = 0.0)
        sqrtI = sqrt(max(I, 0.0))

        # ── Solute log-activities: ln aᵢ = ln10 × log₁₀(γᵢ) + ln(mᵢ) ────────
        for k in eachindex(idx_solutes)
            i       = idx_solutes[k]
            z       = charges[k]
            mᵢ      = max(m[k], ϵ)
            log10_γ = if z == 0
                Kn * I
            else
                -A * z^2 * sqrtI / (1.0 + B * å_vec[k] * sqrtI) + Bdot * I
            end
            out[i] = log10_γ * ln10 + log(mᵢ)
        end

        # ── Water activity via osmotic coefficient ─────────────────────────────
        Σm   = sum(m; init = 0.0)
        Σmz² = 2I                                                       # = Σ mⱼ zⱼ²
        if Σm > 0
            # Charge-squared-weighted mean ion-size parameter
            å_eff = if Σmz² > 0
                sum(
                    m[k] * charges[k]^2 * å_vec[k] for k in eachindex(idx_solutes);
                    init = 0.0,
                ) / Σmz²
            else
                model.å_default
            end
            σ = _hkf_sigma(B * å_eff * sqrtI)
            φ = 1.0 - (A * ln10 / 3) * (Σmz² / Σm) * sqrtI * σ + (Bdot * ln10 / 2) * I
            out[idx_solvent] = -M_w * Σm * φ       # ln a_w
        end

        # ── Gas phase: ideal mixture ───────────────────────────────────────────
        n_gas = sum((_n[i] for i in idx_gas); init = 0.0)
        if !iszero(ustrip(n_gas)) && !isempty(idx_gas)
            for i in idx_gas
                out[i] = log(_n[i] / n_gas)
            end
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
