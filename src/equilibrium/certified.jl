# ── certified.jl ──────────────────────────────────────────────────────────────
#
# An equilibrium that comes with a proof, or says it has none.
#
# The back ends do not agree, and the disagreement is not small. Measured on the
# same problems, with the element balance judged row by row against each row's
# own budget:
#
#   | case                        | interior point | dual Newton |
#   |-----------------------------|----------------|-------------|
#   | calcite in water, 10-40 °C  | 3.0e-2         | 3.0e-12     |
#   | calcite + 1 mmol CO2        | 1.0e-3         | 6.3e-13     |
#   | calcite + 50 mmol CO2       | 1.3e-16        | 8.5e-14     |
#   | pure water                  | 1.3e-16        | 1.3e-16     |
#   | CEM I paste, w/c 0.45, 0.60 | 1.5e-14        | 2.0e-14     |
#   | CEM I paste, w/c 0.30       | 7.2e-16        | 3.3e-10, not certified |
#
# The interior point is wrong in the second digit of the charge balance on a
# calcite solution, because the fraction-to-boundary rule caps its step and the
# residual stops moving; the dual Newton fails to admit a supersaturated phase on
# a low-water cement. Neither is the right default on its own.
#
# What settles it is that the problem is convex whenever the mixing terms are, so
# the KKT conditions are sufficient and `optimality_certificate` DECIDES rather
# than ranks. Given a decision procedure, solving by every available route and
# keeping a proved answer is exact, not heuristic.

"""
    _MAX_RESTARTS

How many times [`equilibrate_certified`](@ref) may restart from its own answer
before giving up. One round is what the measured cases need; the bound exists so
a case that improves by a hair every round cannot loop.
"""
const _MAX_RESTARTS = 3

"""
    _keep_better(eq, cert, eq2, cert2) -> (eq, cert)

Keep the better of two answers: a certificate of optimality beats none, and
otherwise the smaller KKT error wins. A round that buys nothing changes nothing,
which is what lets the restart loop run without ever making the answer worse.

The optimality flag is compared **first**, in both directions. Ranking on the
KKT error alone would let an uncertified point with a smaller stationarity
displace a certified one, and no residual is worth trading a proof for.
"""
function _keep_better(eq, cert, eq2, cert2)
    cert2.optimal == cert.optimal || return cert2.optimal ? (eq2, cert2) : (eq, cert)
    return cert2.stationarity < cert.stationarity ? (eq2, cert2) : (eq, cert)
end

"""
    equilibrate_certified(state; model, ϵ, b, verbose, autostart) -> (state, certificate)

Equilibrium composition together with a proof of its global optimality, obtained
by solving from every registered back end and keeping the answer
[`optimality_certificate`](@ref) proves optimal.

# The starting point is found, not asked for

When no back end certifies from the state as given, an initial approximation is
computed by continuation — [`homotopy_initial_state`](@ref) — and every back end
is run again from it. This is what makes a realistic cement solvable without the
caller knowing anything about the answer: from the cold state of a CEM I paste
(all the mass in the reactants, every product at the `ϵ` floor) no route reaches
the optimum, and with the continuation the same call certifies.

It costs nothing in the ordinary case, because it only runs when nothing else
certified. `autostart = false` declines it, which is what the coupled kinetic
step does: there the caller already supplies the previous instant as a warm
start, and a handful of extra solves inside an implicit ODE step would be paid
at every step.

`certificate.optimal == true` is a **proof**, valid because the Gibbs
minimization is convex when the mixing terms are — ideal mixing and any activity
model whose excess Gibbs energy is convex in the amounts. It is not a proof for a
model that is not, and none of the activity models that ship here have been shown
to violate it; `HKFActivityModel`, `DaviesActivityModel` and the Redlich–Kister
solid solutions are used within their stated ranges.

When no route yields a proof, the answer with the smallest KKT error is returned,
its certificate says so, and a warning names the residual. That is the honest
outcome, and it is not the same thing as a failure: on a low-water cement the
returned composition satisfies the element balance to 1e-15 and has a
supersaturated phase left out, which the certificate reports as
`worst_supersaturation > 0`.

Requires `OptimaSolver` (the dual Newton lives there). Systems without an aqueous
phase, or without `H2O@`, cannot use the dual route; for those, this falls back to
the plain [`equilibrate`](@ref) and returns `nothing` as the certificate.

# Example

```julia
using ChemistryLab, OptimaSolver
eq, cert = equilibrate_certified(state)
cert.optimal          # true — proved globally optimal
cert.balance          # element balance residual
cert.worst_supersaturation   # negative: every absent phase undersaturated
```
"""
function equilibrate_certified(
        state::ChemicalState;
        model::AbstractActivityModel = DiluteSolutionModel(),
        b = nothing,
        ϵ::Float64 = 1.0e-16,
        verbose::Bool = false,
        constraint::EquilibriumConstraint = FixedTP(),
        parameters::Union{Nothing, Base.RefValue} = nothing,
        autostart::Bool = true,
        kwargs...,
    )
    if !_DUAL_AVAILABLE[]
        error(
            "equilibrate_certified needs `OptimaSolver`: the KKT solver and the " *
                "certificate live there. Add `using OptimaSolver` — and " *
                "optionally `using Optimization, OptimizationIpopt` as a second " *
                "starting route, since neither back end certifies every case.",
        )
    end

    if !_dual_applicable(state.system)
        constraint isa FixedTP || throw(
            ArgumentError(
                "a constraint other than `FixedTP` needs the dual route, which " *
                    "requires an aqueous phase and `H2O@` among the species: the " *
                    "prescribed property is an unknown of that solver's system, " *
                    "and the interior-point back ends have nowhere to put it.",
            )
        )
        # No dual route: return the plain answer and say there is no proof.
        eq = equilibrate(state; model = model, ϵ = ϵ, certify = false, kwargs...)
        return (eq, nothing)
    end

    # Duals take the implicit-function route, dispatched on the state's element
    # type rather than tested for. See `_certified_primal_then_derivative`.
    dual_route = _certified_dual_route(
        _amount_number_type(state), state, model, b, ϵ, verbose, constraint,
        parameters, kwargs,
    )
    dual_route === nothing || return dual_route

    des = DualEquilibriumSolver(state.system, model; verbose = verbose)

    # `b` is fixed ONCE, from the state as given. Letting each start define its
    # own would pose a different problem for each: a start that violates the
    # balance — the interior point does, by 3e-6 mol on this class of problem —
    # shifts the component totals by exactly its own infeasibility, and the dual
    # solve then certifies the answer to the shifted problem. Measured, that gave
    # two "certified" compositions 0.2 % apart on dissolved calcium, which on a
    # convex problem with one minimum can only mean two different problems.
    bfix = b === nothing ?
        des.A * Float64[ustrip(us"mol", x) for x in state.n] :
        Float64.(collect(b))

    # Every back end's answer from `from`, and `from` itself — the only start
    # available if they all threw.
    function starts_from(from::ChemicalState, what::AbstractString)
        out = ChemicalState[]
        for f in _SOLVER_FACTORIES
            try
                esolver = EquilibriumSolver(state.system, model, f(); kwargs...)
                push!(out, SciMLBase.solve(esolver, from; ϵ = ϵ, b = bfix))
            catch err
                verbose && @info "$what rejected" backend = f err
            end
        end
        push!(out, from)
        return out
    end

    starts = starts_from(state, "start")

    eq, cert = solve_certified(
        des, starts; b = bfix, ϵ = ϵ, constraint = constraint, parameters = parameters,
    )

    # An automatic initial approximation, computed rather than asked for.
    #
    # Only when nothing above certified, so the common case pays nothing for it.
    # A realistic cement does not converge from the state as given — all the mass
    # in the reactants, every product at the `ϵ` floor — and the caller should
    # not have to know that, nor supply a chemically informed guess.
    # `homotopy_initial_state` walks the solute amount up from a dilute system,
    # which costs a handful of extra solves and needs nothing from the caller.
    if autostart && !cert.optimal
        # Walked under the IDEAL model, deliberately, whatever `model` is: the
        # non-ideal ones do not walk (the a = 0 Debye-Huckel runs away to
        # I = 18 mol/kg, its coefficients falling with I raising solubility
        # raising I). The ideal endpoint is then a good start for `model`,
        # which is what the back-end loop below does with it.
        guess = homotopy_initial_state(state; ϵ = ϵ, verbose = verbose)
        if guess !== nothing
            eq, cert = _keep_better(
                eq, cert,
                solve_certified(
                    des, vcat(starts_from(guess, "start from the continuation"), starts);
                    b = bfix, ϵ = ϵ, constraint = constraint, parameters = parameters,
                )...,
            )
        end

        # Restart from the answer. The continuation ends on a composition that is
        # nearly the equilibrium but not certifiably so, and one more solve from
        # there closes the gap — measured on a CEM I paste under the per-species
        # Debye-Huckel model, stationarity 9.9e-7 (uncertified) becomes 1.5e-16
        # with the worst absent phase 1.4e-5 below saturation. It is the same
        # observation that motivates the continuation, applied once more: a start
        # near the answer is what this problem needs, and the best one available
        # is the answer already in hand.
        #
        # Bounded, and it stops as soon as a round buys nothing, so a genuinely
        # hard case costs a fixed handful of solves rather than looping.
        for _ in 1:_MAX_RESTARTS
            cert.optimal && break
            eq2, cert2 = solve_certified(
                des, starts_from(eq, "restart from the answer"); b = bfix, ϵ = ϵ,
                constraint = constraint, parameters = parameters,
            )
            improved = cert2.optimal || cert2.stationarity < cert.stationarity
            eq, cert = _keep_better(eq, cert, eq2, cert2)
            improved || break
        end
    end

    if !cert.optimal
        # `STRICT_CONVERGENCE[]` is honored here, not only on the interior-point
        # retcode. A caller who sets it is asking that a non-converged solve
        # never pass as a result, and an uncertified answer from this route is
        # exactly that: it can violate the element balance by moles and still
        # come back looking like an ordinary `ChemicalState` — measured, a paste
        # returned with a balance off by 6.7 mol, every hydrate at zero and a
        # table of amounts that reads as a result. A warning is the right default
        # (the answer is still the best one found, and `optimality_certificate`
        # audits it), but under the strict flag it must raise.
        msg = "no route produced a certifiable equilibrium: stationarity " *
            "$(cert.stationarity), element balance $(cert.balance), worst " *
            "supersaturation $(cert.worst_supersaturation)"
        STRICT_CONVERGENCE[] && error(
            msg * ". `ChemistryLab.STRICT_CONVERGENCE[]` is set, so this raises " *
                "rather than returning an answer that is not an equilibrium. " *
                "Audit it with `optimality_certificate`, and see `autostart` for " *
                "the automatic initial approximation."
        )
        @warn msg * "; returning the answer with the smallest KKT error — audit it with `optimality_certificate`" maxlog = 1
    end
    return (eq, cert)
end

"""
    _dual_applicable(system) -> Bool

Whether [`DualEquilibriumSolver`](@ref) can be built for `system`: it needs an
aqueous phase, and `H2O@` among the species, because it parameterizes the interior
variables by the solvent's chemical potential.
"""
function _dual_applicable(system::ChemicalSystem)
    isempty(system.idx_aqueous) && return false
    return haskey(system.dict_species, "H2O@")
end


"""
    _certified_dual_route(state, model, b, ϵ, verbose, constraint, parameters, kwargs)

`nothing` for a real-valued composition; the certified answer with its derivative
attached for one carrying `ForwardDiff.Dual` amounts.

Dispatched on the element type, positionally, so the choice is the type system's
and neither path pays for the other. Two paths are needed for a mathematical
reason, not for want of a generic element type: making the component totals
generic would let duals flow into the solver, and what came back would be the
derivative of the **algorithm** — an active set decided by sign tests, a line
search with branches, an iteration count that varies with the data — rather than
the derivative of the **solution**. The map `b ↦ n*(b)` is smooth only piecewise,
and on each piece the implicit function theorem gives its derivative at the
solution with the active set frozen. That is how `OptimaSolver` computes its own
`Sensitivity`, and how Optima does upstream.
"""
_amount_number_type(state::ChemicalState) = _number_type(eltype(state.n))
_number_type(::Type{<:DynamicQuantities.AbstractQuantity{T}}) where {T} = T
_number_type(::Type{T}) where {T <: Real} = T

_certified_dual_route(::Type{<:Real}, state, model, b, ϵ, verbose, constraint, parameters, kwargs) =
    nothing

function _certified_dual_route(
        ::Type{<:ForwardDiff.Dual}, state, model, b, ϵ, verbose, constraint,
        parameters, kwargs,
    )
    eq_v, cert = equilibrate_certified(
        _primal(state); model = model, ϵ = ϵ, verbose = verbose,
        constraint = constraint, parameters = parameters,
        b = b === nothing ? nothing : _plain.(b), kwargs...,
    )
    nstar = Float64[ustrip(us"mol", x) for x in eq_v.n]
    μ = build_potentials(state.system, model)
    return (_attach_sensitivity(state, nstar, μ, ϵ; b = b), cert)
end
