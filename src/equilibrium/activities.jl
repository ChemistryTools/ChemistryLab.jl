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
