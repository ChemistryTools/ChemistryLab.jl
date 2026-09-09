# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities

"""
    EquilibriumProblem

Definition of a chemical equilibrium problem.

# Fields

  - `b`: conservation vector (elemental abundances).
  - `A`: stoichiometric matrix (conservation matrix).
  - `μ`: chemical potential function `μ(n, p)`.
  - `u0`: initial guess for species amounts.
  - `p`: coefficients for the potential function (default: `nothing`).
  - `lb`: lower bounds for species amounts.
  - `ub`: upper bounds for species amounts.

The problem solves for the species distribution that minimizes the Gibbs energy
subject to mass conservation constraints `A * n = b`.
"""
struct EquilibriumProblem{F <: Function, Tb, TA, Tu, P}
    b::Vector{Tb}
    A::Matrix{TA}
    μ::F
    u0::Vector{Tu}
    p::P
    lb::Vector{Tu}
    ub::Vector{Tu}
end

"""
    _concrete_float(x) -> AbstractArray

Narrow a numeric array to a concrete floating-point element type, returning an
already-concrete one untouched.

`ChemicalSystem` stores its stoichiometry as `Matrix{Real}` whenever integer and
rational coefficients coexist, which a cement's does — `C3AFS0.84H4.32` and its
kind. That is right for the chemistry and wrong for the solver: an abstract
element type boxes every entry and turns `mul!(res, A, x)` into the generic
fallback with a dynamic dispatch per element, on a product evaluated at every
objective and constraint call. SciMLBase says so out loud, warning that "arrays
or dicts to store parameters of different types can hurt performance" as soon as
such an array reaches the problem's parameters, and the warning is correct.

Applies to the default `b = A * u0` as well, which inherits the same abstract
element type from `A`.

Floating point rather than the promoted exact type (`Rational{Int}` here) because
the numeric pipeline downstream is `Float64` throughout — `u0`, the bounds and
`DualEquilibriumSolver.A`, which has always converted — so keeping rationals
would only pay for a rational-to-float conversion at every entry of every
product. The conversion is **lossy for a non-dyadic rational**: `2//5` becomes
`0.4`, which is not equal to it. That is the same rounding the rest of the solver
already applies, and the exact matrix is untouched in `system.SM.A`, where the
stoichiometry belongs.

A caller who passes a concrete array keeps exactly what they passed, identically:
an exact `Matrix{Rational{Int}}`, a `Matrix{Int}`, or a `Matrix{<:Dual}` for
someone differentiating through it.
"""
function _concrete_float(x::AbstractArray)
    # Not restricted to `AbstractArray{<:Number}`: the array that most needs this
    # is the default `b = A * u0`, and `Matrix{Real} * Vector{Float64}` comes out
    # as `Vector{Any}`, whose element type is not a `Number` subtype at all.
    isconcretetype(eltype(x)) && return x
    Tc = float(mapreduce(typeof, promote_type, x; init = Bool))
    return convert(typeof(similar(x, Tc)), x)
end

"""
    EquilibriumProblem(A, μ, u0; b=A*u0, p=nothing, lb=fill(Tu(1e-16), length(u0)), ub=maximum(abs.(A))/minimum(abs.(A[.!iszero.(A)]))*sum(u0)*one.(u0))

Construct an `EquilibriumProblem` with the given stoichiometric matrix `A`, chemical potential function `μ`, and initial guess `u0`.

# Arguments

  - `A`: stoichiometric matrix (conservation matrix).
  - `μ`: chemical potential function `μ(n, p)`.
  - `u0`: initial guess for species amounts.
  - `b`: conservation vector (elemental abundances). Defaults to `A * u0`.
  - `p`: coefficients for the potential function. Defaults to `nothing`.
  - `lb`: lower bounds for species amounts. Defaults to `fill(Tu(1e-16), length(u0))`.
  - `ub`: upper bounds for species amounts. Defaults to `maximum(abs.(A))/minimum(abs.(A[.!iszero.(A)]))*sum(u0)*one.(u0)`.

# Returns

An `EquilibriumProblem` instance.
"""
function EquilibriumProblem(
        A::AbstractMatrix{TA},
        μ::F,
        u0::AbstractVector{Tu};
        b::AbstractVector = A * u0,
        p = nothing,
        lb::AbstractVector = fill(Tu(1.0e-16), length(u0)),
        ub::AbstractVector = maximum(abs.(A)) / minimum(abs.(A[.!iszero.(A)])) * sum(u0) * one.(u0),
    ) where {Tu <: Number, TA <: Number, F <: Function}
    ϵ = 1.0e-16
    lb = max.(lb, ϵ)
    ub = max.(ub, ϵ)
    # Ensure u0 has no zeros or negative values
    u0 = max.(u0, ϵ)
    # `A` arrives with an abstract element type, and the default `b = A * u0`
    # inherits it. See `_concrete_float`.
    Ac = _concrete_float(A)
    bc = _concrete_float(b)
    Tb = eltype(bc)
    return EquilibriumProblem{F, Tb, eltype(Ac), Tu, typeof(p)}(
        Vector{Tb}(bc), Matrix{eltype(Ac)}(Ac), μ, Vector{Tu}(u0), p,
        Vector{Tu}(lb), Vector{Tu}(ub),
    )
end

# Solution transformation methods using multiple dispatch
_solution_transform(::Val{:linear}) = identity
_solution_transform(::Val{:log}) = exp
