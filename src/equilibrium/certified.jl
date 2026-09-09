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
residuals alone would let an uncertified point displace a certified one, and no
residual is worth trading a proof for.

Among uncertified answers the comparison is [`_kkt_error`](@ref) — the worst of
the three residuals — and not the stationarity alone, which is the same ranking
[`solve_certified`](@ref) uses internally. The distinction is not academic: an
answer stationary to 2.4e-3 whose element balance was off by 6.7 mol beat every
candidate the continuation produced, because those were stationary to only 1e-2
while conserving mass. The search computed a usable answer and discarded it.
"""
function _keep_better(eq, cert, eq2, cert2)
    cert2.optimal == cert.optimal || return cert2.optimal ? (eq2, cert2) : (eq, cert)
    return _kkt_error(cert2) < _kkt_error(cert) ? (eq2, cert2) : (eq, cert)
end

"""
    _repair_start(eq, model, bfix, ϵ) -> Union{ChemicalState, Nothing}

Build a starting point that makes the phases the certificate says are missing
actually present. `nothing` when there are none.

This is what turns a diagnosis into a repair. `optimality_certificate` reports a
positive worst supersaturation when a phase sits at the lower bound while the
solution is supersaturated with respect to it — a genuine KKT failure on a convex
problem, so the active set is wrong and the answer is not the answer.
[`saturation_indices`](@ref) says *which* phase, and the obvious move is then to
put it in and solve again.

The failure it exists for is a **phase swap**, which an active-set loop that
admits one phase at a time cannot perform. Measured on a CEM I paste, reported
from one machine while another certified the same source: all 0.02515 mol of
magnesium sat in brucite with `hydrotalcite` absent and supersaturated by 5.58
log units, and admitting the hydrotalcite requires dissolving the brucite
entirely and taking aluminum back from the hydrogarnet in the same step. The
solve was otherwise impeccable — stationarity 1.5e-16, element balance 1.8e-14 —
which is exactly what a converged-onto-the-wrong-active-set answer looks like.

The amount each missing phase is given is what the recipe could make of it,
`min_c b_c / A_cs` over the components it consumes, scaled by `_REPAIR_FRACTION`.
That is a chemical bound, not a guess at the answer: a carbonate in a system with
1e-9 mol of carbon is offered 1e-9 mol and no more. What matters is only that the
phase starts well away from the boundary, since being *at* the boundary is what
the active-set loop cannot recover from. The element balance of the result is not
this function's business — `b` is fixed once by the caller and every start is
projected onto it, so a start that over-spends the budget costs nothing.
"""
function _repair_start(eq::ChemicalState, model, bfix, ϵ::Float64)
    cs = eq.system
    si = saturation_indices(eq, model; ϵ = ϵ)
    n = Float64[ustrip(us"mol", x) for x in eq.n]
    A = cs.SM.A
    # Solid-solution end-members are excluded: their activity is `ln x`, so the
    # index of one sitting at the bound says its mole fraction is small, not that
    # the phase should form. Only pure phases are repaired here.
    ss = Set(cs.idx_ssendmembers)
    missing_phases = [
        i for i in cs.idx_crystal
            if !(i in ss) && n[i] <= 10ϵ && get(si, symbol(cs.species[i]), -Inf) > 1.0e-4
    ]
    isempty(missing_phases) && return nothing

    n2 = copy(n)
    for i in missing_phases
        cap = minimum(
            (
                bfix[c] / A[c, i]
                    for c in eachindex(bfix) if A[c, i] > 0 && bfix[c] > 0
            );
            init = Inf,
        )
        isfinite(cap) && cap > 0 || continue
        n2[i] = max(n2[i], _REPAIR_FRACTION * cap)
    end
    n2 == n && return nothing
    return ChemicalState(
        cs, n2 .* u"mol"; T = temperature(eq), P = pressure(eq),
    )
end

"""
    _REPAIR_FRACTION

What fraction of the amount the recipe could make of a missing phase
[`_repair_start`](@ref) puts in. Far enough off the boundary for the active-set
loop to work with, small enough not to pretend it knows the answer.
"""
const _REPAIR_FRACTION = 0.1

"""
    _repair_round(eq, cert, model, bfix, ϵ, solve_from, verbose)
        -> (eq, cert, improved)

One round of "the certificate names a missing phase, so put it in and solve
again". Returns the better of the two answers and whether the round bought
anything, so the caller can stop as soon as it does not.

`solve_from` is the search, injected: given a starting composition it returns
`(state, certificate)`. Passing it in rather than closing over the caller's
solver is what makes this round testable on its own — the situation it exists
for, a back end that converges onto the wrong active set, needs a system of some
135 species to arise, while the round's logic needs eight.
"""
function _repair_round(eq, cert, model, bfix, ϵ::Float64, solve_from, verbose::Bool)
    fixed = _repair_start(eq, model, bfix, ϵ)
    fixed === nothing && return (eq, cert, false)
    verbose && @info "repairing a missing phase" worst_si = cert.worst_supersaturation
    eq2, cert2 = solve_from(fixed)
    improved = cert2.optimal || _kkt_error(cert2) < _kkt_error(cert)
    kept, kept_cert = _keep_better(eq, cert, eq2, cert2)
    return (kept, kept_cert, improved)
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
    #
    # `STRICT_CONVERGENCE[]` is cleared for the duration and restored after: what
    # this computes is a STARTING POINT, not a result. Left set, a back end that
    # reports `MaxIters` raises, the `catch` below swallows it, and the search
    # silently loses that candidate — so a caller asking for strict results gets
    # a *worse* search than a caller who did not. Measured on a CEM I paste where
    # the interior point ends on `MaxIters`: with the flag set the route returned
    # an element balance of 27.6 mol, and with it clear the very same call
    # returned 1.8e-14. The result is still judged strictly, at the end of this
    # function, which is where the flag belongs.
    function starts_from(from::ChemicalState, what::AbstractString)
        out = ChemicalState[]
        strict = STRICT_CONVERGENCE[]
        STRICT_CONVERGENCE[] = false
        try
            _exploring_starts() do
                for f in _SOLVER_FACTORIES
                    try
                        esolver = EquilibriumSolver(state.system, model, f(); kwargs...)
                        push!(out, SciMLBase.solve(esolver, from; ϵ = ϵ, b = bfix))
                    catch err
                        verbose && @info "$what rejected" backend = f err
                    end
                end
            end
        finally
            STRICT_CONVERGENCE[] = strict
        end
        push!(out, from)
        return out
    end

    # Every candidate the search tries is a candidate, and a candidate that does
    # not converge is what the search exists for. Its diagnostics stay quiet; the
    # verdict on the answer is pronounced once, below, on the certificate.
    search(starts) = _exploring_starts() do
        solve_certified(
            des, starts; b = bfix, ϵ = ϵ,
            constraint = constraint, parameters = parameters,
        )
    end

    starts = starts_from(state, "start")

    eq, cert = search(starts)

    # An automatic initial approximation, computed rather than asked for.
    #
    # Only when nothing above certified, so the common case pays nothing for it.
    # A realistic cement does not converge from the state as given — all the mass
    # in the reactants, every product at the `ϵ` floor — and the caller should
    # not have to know that, nor supply a chemically informed guess.
    # `homotopy_initial_state` walks the solute amount up from a dilute system,
    # which costs a handful of extra solves and needs nothing from the caller.
    # What the automatic start did, in words, for the diagnostic below. When a
    # solve fails on one machine and not another, the first thing anyone needs to
    # know is whether the continuation ran at all and whether it helped — and
    # asking for that should not require a second run with `verbose = true`.
    note = autostart ? "not reached (the first route certified)" : "declined (autostart = false)"
    if autostart && !cert.optimal
        # Walked under the IDEAL model, deliberately, whatever `model` is: the
        # non-ideal ones do not walk (the a = 0 Debye-Huckel runs away to
        # I = 18 mol/kg, its coefficients falling with I raising solubility
        # raising I). The ideal endpoint is then a good start for `model`,
        # which is what the back-end loop below does with it.
        guess = homotopy_initial_state(state; ϵ = ϵ, verbose = verbose)
        if guess === nothing
            note = "the continuation produced no usable start: every rung was " *
                "refused, or the system has no aqueous solvent to walk"
        else
            before = _kkt_error(cert)
            eq, cert = _keep_better(
                eq, cert,
                search(vcat(starts_from(guess, "start from the continuation"), starts))...,
            )
            note = cert.optimal ?
                "the continuation certified it" :
                (
                    _kkt_error(cert) < before ?
                    "the continuation improved the KKT error from $before to " *
                    "$(_kkt_error(cert)) without certifying" :
                    "the continuation ran and its answer was no better than " *
                    "$before, so it was discarded"
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
            eq2, cert2 = search(starts_from(eq, "restart from the answer"))
            improved = cert2.optimal || _kkt_error(cert2) < _kkt_error(cert)
            eq, cert = _keep_better(eq, cert, eq2, cert2)
            improved || break
        end

        # Act on what the certificate says. A positive worst supersaturation
        # names a phase that should be present and is not, which no amount of
        # restarting from the same active set will fix: see `_repair_start`.
        # The accumulated starts go back in with the repaired composition.
        # Measured on the paste, a solve from the repair start alone reaches the
        # right assemblage -- 74.19 cm3, all the magnesium back in the
        # hydrotalcite -- and still fails its certificate on an unrelated trace
        # component: the recipe's 1e-9 mol of carbon is lost, leaving an element
        # balance of exactly 1e-9 against a tolerance of 1e-10. Ranked on the
        # worst residual, that answer loses to the very point it was meant to
        # replace. Handing the search the repaired composition *and* the
        # candidates it already had keeps the chemistry of the one and the trace
        # components of the others.
        repair_search(f) = search(vcat(starts_from(f, "repair start"), starts))
        for _ in 1:_MAX_RESTARTS
            cert.optimal && break
            eq, cert, improved = _repair_round(
                eq, cert, model, bfix, ϵ, repair_search, verbose,
            )
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
            "supersaturation $(cert.worst_supersaturation). Automatic initial " *
            "approximation: $note"
        STRICT_CONVERGENCE[] && error(
            msg * ". `ChemistryLab.STRICT_CONVERGENCE[]` is set, so this raises " *
                "rather than returning an answer that is not an equilibrium. " *
                "Audit it with `optimality_certificate`; " *
                "`homotopy_initial_state(state; verbose = true)` reports each rung."
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
