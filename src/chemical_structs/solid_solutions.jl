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
state and `SC_SSENDMEMBER` class) mixing according to `model`.

# Construction

Use the keyword constructor:
```julia
SolidSolutionPhase(name, end_members; model = IdealSolidSolutionModel())
```

Validation at construction time:
- All end-members must have `aggregate_state == AS_CRYSTAL` and `class == SC_SSENDMEMBER`.
- [`RedlichKisterModel`](@ref) requires exactly 2 end-members.

# Example

```jldoctest
julia> em1 = Species("Ca2SiO4"; aggregate_state=AS_CRYSTAL, class=SC_SSENDMEMBER);

julia> em2 = Species("Ca3Si2O7"; aggregate_state=AS_CRYSTAL, class=SC_SSENDMEMBER);

julia> ss = SolidSolutionPhase("CSH", [em1, em2])
SolidSolutionPhase{Species{Int64}, IdealSolidSolutionModel}
  name: CSH
  end-members (2): Ca2SiO4, Ca3Si2O7
  model: IdealSolidSolutionModel
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
        class(sp) == SC_SSENDMEMBER ||
            error(
            "SolidSolutionPhase: end-member \"$(symbol(sp))\" must have " *
                "class = SC_SSENDMEMBER (got $(class(sp)))",
        )
    end
    if model isa RedlichKisterModel
        length(end_members) == 2 ||
            error(
            "RedlichKisterModel requires exactly 2 end-members, " *
                "got $(length(end_members))",
        )
    end
    T = eltype(end_members)
    return SolidSolutionPhase{T, typeof(model)}(
        String(name), collect(T, end_members), model
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
