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

    # Resolve indices once at construction time — not at evaluation time
    idx_solvent = only(cs.idx_solvent)          # unique solvent index
    idx_solutes = cs.idx_solutes                # aqueous solute indices
    idx_crystal = cs.idx_crystal                # pure solid indices
    idx_gas     = cs.idx_gas                    # gas species indices

    # Molar concentration of pure solvent c°_solvent = 1/M_solvent [mol/L]
    # Stripped of units at construction — used as a constant in the closure
    M_solvent = ustrip(us"kg/mol", cs.species[idx_solvent][:M])  # kg/mol
    ρ_solvent = 1.0                              # kg/L for water — standard reference
    c_solvent = ρ_solvent / M_solvent            # mol/L — ~55.5 for water

    # Precompute correction for solutes: ln(c°/c_solvent) = -ln(c_solvent)
    # Henry convention: ln a_i = ln(n_i/n_solvent) + ln(c_solvent)
    ln_c_solvent = log(c_solvent)

    function lna(n::AbstractVector, p)
        ϵ   = p.ϵ                               # floor to avoid log(0)
        _n  = max.(n, ϵ)                        # regularized mole vector

        lna = zeros(eltype(_n), length(_n))     # supports Dual numbers for autodiff

        # Solvent — Raoult: ln a = ln(n_solvent / n_aqueous_total)
        n_aqueous = _n[idx_solvent] + sum(_n[idx_solutes])
        if !iszero(n_aqueous)
            lna[idx_solvent] = log(_n[idx_solvent] / n_aqueous)
        end

        # Solutes — Henry: ln a_i = ln(n_i / n_solvent) + ln(c_solvent)
        # Equivalent to ln(c_i [mol/L]) with c° = 1 mol/L reference
        if !iszero(_n[idx_solvent])
            for i in idx_solutes
                lna[i] = log(_n[i] / _n[idx_solvent]) + ln_c_solvent
            end
        end

        # Crystals — pure solid: a = 1, ln a = 0 (already zero)

        # Gas — ideal mixture: ln a_i = ln(x_i) = ln(n_i / n_gas_total)
        n_gas = sum(_n[idx_gas])
        if !iszero(n_gas)
            for i in idx_gas
                lna[i] = log(_n[i] / n_gas)
            end
        end

        return lna
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
