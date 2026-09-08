# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities
using ForwardDiff
using OrderedCollections

# ── Solid solution activity models ────────────────────────────────────────────

"""
    abstract type AbstractSolidSolutionModel end

Base type for activity models within a solid solution phase.
Each concrete subtype implements `_excess_ln_gamma(model, k, x, T)`.
"""
abstract type AbstractSolidSolutionModel end

"""
    struct IdealSolidSolutionModel <: AbstractSolidSolutionModel

Ideal (Temkin) solid solution mixing: `ln γᵢ = 0`, hence `ln aᵢ = ln xᵢ`.

Valid for any number of end-members.

# Example

```jldoctest
julia> m = IdealSolidSolutionModel()
IdealSolidSolutionModel()
```
"""
struct IdealSolidSolutionModel <: AbstractSolidSolutionModel end

"""
    struct RegularSolutionModel{T<:Real} <: AbstractSolidSolutionModel

Symmetric regular (multi-component Margules) solid solution.

The gap this fills is arity. [`RedlichKisterModel`](@ref) is the general binary
form and is restricted to two end-members, while [`IdealSolidSolutionModel`](@ref)
takes any number but no interaction at all. A C-S-H with six end-members, or the
CNASH and ECSH families of CEMDATA18, had no non-ideal option here.

Excess Gibbs energy, with one symmetric interaction parameter per pair:

```
G^ex = Σ_{i<j} W_ij x_i x_j
ln γ_k = (1/RT) [ Σ_{j≠k} W_kj x_j  −  Σ_{i<j} W_ij x_i x_j ]
```

`W` is in J/mol, symmetric, with a zero diagonal; only the off-diagonal entries
are read and `W[i,j]` is used for the pair `(i,j)`. `W_ij > 0` is repulsive —
it favors unmixing — and `W_ij = 0` recovers ideal mixing.

For two end-members this is exactly `RedlichKisterModel(a0 = W₁₂)`, which the
test suite checks; the point of it is `n > 2`.

AD-compatible: all computations propagate `ForwardDiff.Dual` numbers.

# Examples

```jldoctest
julia> m = RegularSolutionModel([0.0 4000.0; 4000.0 0.0]);

julia> m.W[1, 2]
4000.0
```

# References

  - Guggenheim, E.A. (1937). *Trans. Faraday Soc.* **33**, 151–159.
"""
struct RegularSolutionModel{T <: Real} <: AbstractSolidSolutionModel
    W::Matrix{T}
end

"""
    RegularSolutionModel(W::AbstractMatrix) -> RegularSolutionModel

Construct a [`RegularSolutionModel`](@ref) from a symmetric matrix of
interaction parameters in J/mol. Raises if `W` is not square or not symmetric;
the diagonal is ignored.
"""
function RegularSolutionModel(W::AbstractMatrix)
    n, m = size(W)
    n == m || throw(ArgumentError("RegularSolutionModel: W must be square, got $(size(W))"))
    for i in 1:n, j in (i + 1):n
        isapprox(W[i, j], W[j, i]; rtol = 1.0e-12) || throw(
            ArgumentError(
                "RegularSolutionModel: W must be symmetric; W[$i,$j] = $(W[i, j]) " *
                    "but W[$j,$i] = $(W[j, i]).",
            )
        )
    end
    Wf = float.(collect(W))
    return RegularSolutionModel{eltype(Wf)}(Wf)
end

"""
    struct RedlichKisterModel{T<:Real} <: AbstractSolidSolutionModel

Binary Redlich-Kister (asymmetric Margules) model for non-ideal solid solutions.
Requires **exactly 2 end-members** per solid solution phase.

Parameters `a0`, `a1`, `a2` are in J/mol and divided by RT inside the activity closure.

Activity coefficients (Guggenheim / ThermoCalc convention):

```
ln γ₁ = (x₂²/RT)[a₀ + a₁(3x₁ − x₂) + a₂(x₁ − x₂)(5x₁ − x₂)]
ln γ₂ = (x₁²/RT)[a₀ − a₁(3x₂ − x₁) + a₂(x₂ − x₁)(5x₂ − x₁)]
```

AD-compatible: all computations propagate `ForwardDiff.Dual` numbers.

# Examples

```jldoctest
julia> m = RedlichKisterModel(a0 = 4000.0, a1 = 500.0)
RedlichKisterModel{Float64}(4000.0, 500.0, 0.0)

julia> m.a0
4000.0
```
"""
struct RedlichKisterModel{T <: Real} <: AbstractSolidSolutionModel
    a0::T   # J/mol — symmetric interaction parameter
    a1::T   # J/mol — asymmetry
    a2::T   # J/mol — higher-order correction
end

"""
    RedlichKisterModel(; a0=0.0, a1=0.0, a2=0.0)

Keyword constructor for [`RedlichKisterModel`](@ref). Parameters are promoted to a
common type.
"""
function RedlichKisterModel(; a0 = 0.0, a1 = 0.0, a2 = 0.0)
    vals = promote(a0, a1, a2)
    return RedlichKisterModel{eltype(vals)}(vals...)
end

# ── Solid solution phase ───────────────────────────────────────────────────────

"""
    abstract type AbstractSolidSolutionPhase end

Base type for solid solution phases. Concrete subtypes group a set of end-member
[`AbstractSpecies`](@ref) and associate them with an [`AbstractSolidSolutionModel`](@ref).
"""
abstract type AbstractSolidSolutionPhase end

"""
    struct SolidSolutionPhase{T<:AbstractSpecies, M<:AbstractSolidSolutionModel}
            <: AbstractSolidSolutionPhase

A solid-solution phase consisting of `end_members` (species with `AS_CRYSTAL` aggregate
state) mixing according to `model`.

End-members are automatically requalified to `SC_SSENDMEMBER` at construction time,
so database species with `SC_COMPONENT` can be passed directly.

# Construction

Use the keyword constructor:
```julia
SolidSolutionPhase(name, end_members; model = IdealSolidSolutionModel())
```

Validation at construction time:
- All end-members must have `aggregate_state == AS_CRYSTAL`.
- [`RedlichKisterModel`](@ref) requires exactly 2 end-members.

# Example

```jldoctest
julia> em1 = Species("Ca2SiO4"; aggregate_state=AS_CRYSTAL, class=SC_COMPONENT);

julia> em2 = Species("Ca3Si2O7"; aggregate_state=AS_CRYSTAL, class=SC_COMPONENT);

julia> ss = SolidSolutionPhase("CSH", [em1, em2])
SolidSolutionPhase{Species{Int64}, IdealSolidSolutionModel}
  name: CSH
  end-members (2): Ca2SiO4, Ca3Si2O7
  model: IdealSolidSolutionModel

julia> class(end_members(ss)[1])
SC_SSENDMEMBER::Class = 5
```
"""
struct SolidSolutionPhase{T <: AbstractSpecies, M <: AbstractSolidSolutionModel} <:
    AbstractSolidSolutionPhase
    name::String
    end_members::Vector{T}
    model::M
end

"""
    SolidSolutionPhase(name, end_members; model=IdealSolidSolutionModel())

Construct and validate a [`SolidSolutionPhase`](@ref).

End-members whose `class` is not already `SC_SSENDMEMBER` are automatically
requalified via [`with_class`](@ref). Passing database species with
`SC_COMPONENT` therefore works directly, without a prior call to `with_class`.
"""
function SolidSolutionPhase(
        name::AbstractString,
        end_members::AbstractVector{<:AbstractSpecies};
        model::AbstractSolidSolutionModel = IdealSolidSolutionModel(),
    )
    for sp in end_members
        aggregate_state(sp) == AS_CRYSTAL ||
            error(
            "SolidSolutionPhase: end-member \"$(symbol(sp))\" must have " *
                "aggregate_state = AS_CRYSTAL (got $(aggregate_state(sp)))",
        )
    end
    if model isa RedlichKisterModel
        length(end_members) == 2 ||
            error(
            "RedlichKisterModel requires exactly 2 end-members, " *
                "got $(length(end_members))",
        )
    end
    qualified = [
        class(sp) == SC_SSENDMEMBER ? sp : with_class(sp, SC_SSENDMEMBER)
            for sp in end_members
    ]
    T = eltype(qualified)
    return SolidSolutionPhase{T, typeof(model)}(
        String(name), collect(T, qualified), model
    )
end

# ── Accessors ─────────────────────────────────────────────────────────────────

"""
    name(ss::SolidSolutionPhase) -> String

Return the name of the solid solution phase.
"""
name(ss::SolidSolutionPhase) = ss.name

"""
    end_members(ss::SolidSolutionPhase) -> Vector{<:AbstractSpecies}

Return the end-member species of the solid solution phase.
"""
end_members(ss::SolidSolutionPhase) = ss.end_members

"""
    model(ss::SolidSolutionPhase) -> AbstractSolidSolutionModel

Return the activity model of the solid solution phase.
"""
model(ss::SolidSolutionPhase) = ss.model

# ── Display ───────────────────────────────────────────────────────────────────

function Base.show(io::IO, ss::SolidSolutionPhase{T, M}) where {T, M}
    em_names = join(symbol.(ss.end_members), ", ")
    println(io, "SolidSolutionPhase{$T, $M}")
    println(io, "  name: $(ss.name)")
    println(io, "  end-members ($(length(ss.end_members))): $em_names")
    return print(io, "  model: $M")
end
