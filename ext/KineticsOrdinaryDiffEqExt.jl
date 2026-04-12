"""
    KineticsOrdinaryDiffEqExt

Extension activated when `OrdinaryDiffEq` is loaded. Provides the concrete
`ChemistryLab.integrate(::KineticsProblem, ::KineticsSolver; ...)` implementation
using `OrdinaryDiffEq.ODEProblem` and registers `Rodas5P()` as the default solver.

# Usage

```julia
using ChemistryLab, OrdinaryDiffEq
sol = integrate(kp, KineticsSolver(; ode_solver=Rodas5P()))
# or shortcut (uses default Rodas5P):
sol = integrate(kp)
```
"""
module KineticsOrdinaryDiffEqExt

using OrdinaryDiffEq
import ChemistryLab:
    integrate,
    KineticsProblem,
    KineticsSolver,
    build_kinetics_ode,
    build_u0,
    build_kinetics_params,
    extend_u0,
    extend_ode!,
    n_extra_states,
    _DEFAULT_KINETICS_SOLVER_FACTORY,
    AbstractCalorimeter,
    _total_enthalpy

# ── Concrete integrate implementation ────────────────────────────────────────

"""
    integrate(kp::KineticsProblem, ks::KineticsSolver;
              calorimeter=nothing) -> ODESolution

Integrate the kinetics ODE using `OrdinaryDiffEq`.

The ODE right-hand-side is built by [`build_kinetics_ode`](@ref) and
optionally extended with calorimeter equations via [`extend_ode!`](@ref).

## Algorithm

1. Build `u0` (kinetic mineral moles, optionally extended by calorimeter).
2. Build the immutable parameter `NamedTuple` `p` from the initial state.
3. Wrap into `ODEProblem` and call `OrdinaryDiffEq.solve`.

## Default tolerances

`reltol = 1e-8`, `abstol = 1e-10` — suitable for mineral dissolution kinetics
where mole amounts span several orders of magnitude.

## Returns

A `SciMLBase.ODESolution`. Access mineral moles at time `t` via `sol(t)[i]`
for the i-th kinetic mineral.

# Examples

```julia
using ChemistryLab, OrdinaryDiffEq

ks  = KineticsSolver(; ode_solver=Rodas5P(), reltol=1e-8, abstol=1e-10)
sol = integrate(kp, ks)

# With isothermal calorimetry
using ChemistryLab: IsothermalCalorimeter
cal = IsothermalCalorimeter(298.15)
sol = integrate(kp, ks; calorimeter=cal)
t, Q = cumulative_heat(sol, cal)
```
"""
function integrate(
        kp::KineticsProblem,
        ks::KineticsSolver;
        calorimeter::Union{Nothing, AbstractCalorimeter} = nothing,
    )
    # ── Build base ODE ingredients ─────────────────────────────────────────
    f_base! = build_kinetics_ode(kp)
    u0_base = build_u0(kp)
    p = build_kinetics_params(kp)
    n_kin = length(u0_base)

    # ── Optionally extend state vector for calorimeter ─────────────────────
    u0 = isnothing(calorimeter) ? u0_base : extend_u0(u0_base, calorimeter)

    # Wrap f_base! to also call extend_ode! when a calorimeter is present
    f! = if isnothing(calorimeter)
        f_base!
    else
        cal = calorimeter
        function (du, u, p_, t_)
            f_base!(du, u, p_, t_)
            extend_ode!(du, u, p_, n_kin, cal)
            return nothing
        end
    end

    # ── DiscreteCallback: track total system enthalpy at each accepted step ──
    # Fires at t=0 (initialize) and after every accepted step.
    # Stores H(t) = Σᵢ nᵢ(t)·hᵢ(T) in p.saved_H so that cumulative_heat can
    # return Q(t) = H(0) - H(t), capturing both kinetic and equilibrium heat.
    function _enthalpy_affect!(integrator)
        p_ = integrator.p
        any(!isnothing, p_.h_fns) || return   # no enthalpy data → skip
        # Refresh mineral moles in n_full from the current (accepted) ODE state
        for (j, idx) in enumerate(p_.idx_kin)
            p_.n_full[idx] = max(integrator.u[j], p_.ϵ)
        end
        H = _total_enthalpy(p_.n_full, p_.h_fns, p_.T)
        push!(p_.saved_H, H)
        push!(p_.saved_t, integrator.t)
        return nothing
    end

    enthalpy_cb = DiscreteCallback(
        Returns(true),
        _enthalpy_affect!;
        initialize = (c, u, t, integrator) -> _enthalpy_affect!(integrator),
        save_positions = (false, false),
    )

    # ── Merge default tolerances with user kwargs ──────────────────────────
    defaults = (reltol = 1.0e-8, abstol = 1.0e-10)
    merged = merge(defaults, ks.kwargs)

    solver = isnothing(ks.ode_solver) ? Rodas5P() : ks.ode_solver

    prob = ODEProblem(f!, u0, kp.tspan, p)
    return solve(prob, solver; callback = enthalpy_cb, merged...)
end

# ── __init__: register default solver ────────────────────────────────────────

function __init__()
    return _DEFAULT_KINETICS_SOLVER_FACTORY[] =
        () -> KineticsSolver(; ode_solver = Rodas5P())
end

end  # module KineticsOrdinaryDiffEqExt
