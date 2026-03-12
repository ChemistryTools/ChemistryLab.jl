# ── Abstract mixing model ──────────────────────────────────────────────────────

"""
    abstract type AbstractSolidSolutionMixingModel end

Base type for solid solution mixing models.

Each concrete subtype must implement:
```julia
lnγ_ss(model, x::AbstractVector, RT::Real) -> Vector
```
returning the vector of log activity coefficients `ln γᵢ` for mole fractions
`x` and thermal energy `RT` [J/mol].
"""
abstract type AbstractSolidSolutionMixingModel end

# ── Ideal mixing ───────────────────────────────────────────────────────────────

"""
    struct IdealSolidSolutionMixingModel <: AbstractSolidSolutionMixingModel

Ideal solid solution: `ln γᵢ = 0` for all end-members.
Activity is purely entropic: `ln aᵢ = ln xᵢ`.
"""
struct IdealSolidSolutionMixingModel <: AbstractSolidSolutionMixingModel end

"""
    lnγ_ss(::IdealSolidSolutionMixingModel, x, RT) -> Vector

Ideal activity coefficients: zero vector.
"""
lnγ_ss(::IdealSolidSolutionMixingModel, x::AbstractVector, RT::Real) =
    zeros(eltype(x), length(x))

# ── Symmetric Margules (regular solution) ─────────────────────────────────────

"""
    struct RegularSolidSolutionMixingModel <: AbstractSolidSolutionMixingModel

Symmetric Margules (regular solution) model for a solid solution with n end-members.

## Definition

The excess Gibbs energy per mole of phase is:

```math
G^{xs} = \\sum_{i<j} W_{ij}\\, x_i\\, x_j
```

Activity coefficients are the exact derivative of ``n G^{xs}`` with respect to
``n_i`` at constant ``n_{total}``:

```math
RT \\ln\\gamma_i =
    \\sum_{j \\neq i} W_{ij}\\, x_j (1 - x_i)
  - \\sum_{\\substack{j < k \\\\ j \\neq i,\\, k \\neq i}} W_{jk}\\, x_j\\, x_k
```

For a binary (1, 2) this simplifies to: ``RT \\ln\\gamma_1 = W_{12}\\,x_2^2``.

## Field

  - `W`: ``n \\times n`` symmetric matrix of Margules interaction parameters [J/mol].
    Zero diagonal (`W[i,i] = 0`).

## Constructor

Validates that `W` is square, symmetric, and has a zero diagonal.
"""
struct RegularSolidSolutionMixingModel <: AbstractSolidSolutionMixingModel
    W::Matrix{Float64}  # [J/mol], symmetric, zero diagonal

    function RegularSolidSolutionMixingModel(W::AbstractMatrix{<:Real})
        size(W, 1) == size(W, 2) ||
            throw(ArgumentError("W must be a square matrix"))
        issymmetric(W) ||
            throw(ArgumentError("W must be symmetric (Wᵢⱼ = Wⱼᵢ)"))
        all(iszero(W[i, i]) for i in axes(W, 1)) ||
            throw(ArgumentError("Diagonal of W must be zero (Wᵢᵢ = 0)"))
        return new(Float64.(W))
    end
end

"""
    lnγ_ss(model::RegularSolidSolutionMixingModel, x, RT) -> Vector

Symmetric Margules activity coefficients for mole fractions `x` and thermal
energy `RT` [J/mol].

## Formula

```math
RT \\ln\\gamma_i =
    \\sum_{j \\neq i} W_{ij}\\, x_j (1 - x_i)
  - \\sum_{\\substack{j < k \\\\ j \\neq i,\\, k \\neq i}} W_{jk}\\, x_j\\, x_k
```
"""
function lnγ_ss(model::RegularSolidSolutionMixingModel, x::AbstractVector, RT::Real)
    W   = model.W
    n   = length(x)
    lnγ = zeros(eltype(x), n)
    for i in 1:n
        s = zero(eltype(x))
        for j in 1:n
            j == i && continue
            s += W[i, j] * x[j] * (1 - x[i])
        end
        for j in 1:n
            j == i && continue
            for k in (j + 1):n
                k == i && continue
                s -= W[j, k] * x[j] * x[k]
            end
        end
        lnγ[i] = s / RT
    end
    return lnγ
end

# ── Solid solution descriptor ─────────────────────────────────────────────────

"""
    struct SolidSolution

Groups a set of crystalline end-members (`AS_CRYSTAL`) from a `ChemicalSystem`
and associates them with a mixing model.

## Fields

  - `end_member_symbols`: species symbols of the end-members (must match
    `symbol(sp)` for species in the `ChemicalSystem`).
  - `mixing_model`: mixing model (`IdealSolidSolutionMixingModel`,
    `RegularSolidSolutionMixingModel`, …).

## Constructor

Validates:
- The list must not be empty.
- For `RegularSolidSolutionMixingModel`, `size(W)` must match the number
  of end-members.

Validation that symbols exist in a `ChemicalSystem` and have
`aggregate_state == AS_CRYSTAL` is deferred to `activity_model`.

## Examples

```julia
ss = SolidSolution(["Calcite", "Magnesite"], IdealSolidSolutionMixingModel())

W  = [0.0 6700.0; 6700.0 0.0]
ss = SolidSolution(["Calcite", "Magnesite"], RegularSolidSolutionMixingModel(W))
```
"""
struct SolidSolution
    end_member_symbols::Vector{String}
    mixing_model::AbstractSolidSolutionMixingModel

    function SolidSolution(
        symbols::AbstractVector{<:AbstractString},
        mixing_model::AbstractSolidSolutionMixingModel,
    )
        isempty(symbols) &&
            throw(ArgumentError("SolidSolution requires at least one end-member"))
        if mixing_model isa RegularSolidSolutionMixingModel
            n = length(symbols)
            size(mixing_model.W, 1) == n ||
                throw(
                    DimensionMismatch(
                        "W size $(size(mixing_model.W)) does not match $n end-member(s)",
                    ),
                )
        end
        return new(collect(String, symbols), mixing_model)
    end
end

# ── Composite activity model ──────────────────────────────────────────────────

"""
    struct SolidSolutionActivityModel <: AbstractActivityModel

Decorator that wraps a base activity model and adds solid solution mixing
contributions for designated groups of crystalline end-members.

The base model handles aqueous species, gas species, and pure crystals
(for which `ln a = 0`). This decorator overrides `ln a` for each end-member
in every `SolidSolution` with `ln xᵢ + ln γᵢ`.

## Fields

  - `base`: base activity model (`DiluteSolutionModel`, `HKFActivityModel`, …).
  - `solid_solutions`: vector of `SolidSolution` descriptors.

## Example

```julia
ss = SolidSolution(["Calcite", "Magnesite"], IdealSolidSolutionMixingModel())
model = SolidSolutionActivityModel(DiluteSolutionModel(), [ss])
solver = EquilibriumSolver(cs, model, IpoptOptimizer())
```
"""
struct SolidSolutionActivityModel <: AbstractActivityModel
    base::AbstractActivityModel
    solid_solutions::Vector{SolidSolution}
end

# ── activity_model factory ────────────────────────────────────────────────────

"""
    activity_model(cs::ChemicalSystem, model::SolidSolutionActivityModel) -> Function

Return a closure `lna(n, p) -> Vector` of log-activities.

**At construction time (once):**
- Builds the base model closure.
- Validates that every end-member symbol exists in `cs`, has
  `aggregate_state == AS_CRYSTAL`, and does not belong to more than one group.
- Precomputes integer indices of end-members in `cs.species`.

**At runtime (called every optimizer iteration):**
- Starts from the base closure output (aqueous/gas correct, crystals = 0).
- For each SS group: computes mole fractions `xᵢ`, calls `lnγ_ss`, then
  writes `ln xᵢ + ln γᵢ` at the corresponding indices.
- If the phase is entirely absent (`n_SS ≈ 0`), the group is skipped
  (crystals keep `ln a = 0`).

`p` must contain:
- `ϵ`: regularization floor (avoids `log(0)`).
- `T`: temperature in K (provided by `_build_params`).
"""
function activity_model(cs::ChemicalSystem, model::SolidSolutionActivityModel)

    lna_base = activity_model(cs, model.base)

    # -- Validate and precompute end-member indices (done once at construction) --
    all_ss_indices = Set{Int}()
    ss_index_groups = Vector{Vector{Int}}(undef, length(model.solid_solutions))

    for (k, ss) in enumerate(model.solid_solutions)
        idxs = Int[]
        for sym in ss.end_member_symbols
            haskey(cs.dict_species, sym) ||
                throw(
                    ArgumentError(
                        "Solid solution end-member \"$sym\" not found in ChemicalSystem",
                    ),
                )
            sp = cs.dict_species[sym]
            aggregate_state(sp) == AS_CRYSTAL ||
                throw(
                    ArgumentError(
                        "End-member \"$sym\" must have aggregate_state == AS_CRYSTAL " *
                        "(got: $(aggregate_state(sp)))",
                    ),
                )
            i = findfirst(s -> symbol(s) == sym, cs.species)
            i in all_ss_indices &&
                throw(
                    ArgumentError(
                        "End-member \"$sym\" belongs to more than one SolidSolution",
                    ),
                )
            push!(all_ss_indices, i)
            push!(idxs, i)
        end
        ss_index_groups[k] = idxs
    end

    mixing_models = [ss.mixing_model for ss in model.solid_solutions]

    # -- Closure (hot path — called every optimizer iteration) --
    function lna(n::AbstractVector, p)
        out = lna_base(n, p)  # aqueous/gas correct; crystals = 0.0
        ϵ  = p.ϵ
        RT = ustrip(us"J/mol/K", Constants.R) * p.T  # [J/mol]

        for (k, idxs) in enumerate(ss_index_groups)
            n_SS = sum(max(n[i], ϵ) for i in idxs)
            # Phase entirely absent: keep ln a = 0 from base model
            n_SS <= length(idxs) * ϵ * 10 && continue

            x   = [max(n[i], ϵ) / n_SS for i in idxs]
            lnγ = lnγ_ss(mixing_models[k], x, RT)

            for (j, i) in enumerate(idxs)
                out[i] = log(x[j]) + lnγ[j]
            end
        end

        return out
    end

    return lna
end
