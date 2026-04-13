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
    SemiAdiabaticCalorimeter,
    _total_enthalpy

# ── Concrete integrate implementation ────────────────────────────────────────

"""
    integrate(kp::KineticsProblem, ks::KineticsSolver;
              calorimeter=nothing) -> ODESolution

Integrate the kinetics ODE using `OrdinaryDiffEq`.

## Algorithm

1. Build `u0` (kinetic mineral moles, optionally extended by calorimeter).
2. Build the immutable parameter `NamedTuple` `p` from the initial state.
3. Wrap into `ODEProblem` and call `OrdinaryDiffEq.solve`.

Default tolerances: `reltol = 1e-8`, `abstol = 1e-10`.

When a [`SemiAdiabaticCalorimeter`](@ref) is used, a warning is emitted for
any species that lack Cp° data (their contribution to `Cp_total` is zero).

# Examples

```julia
using ChemistryLab, OrdinaryDiffEq

ks  = KineticsSolver(; ode_solver=Rodas5P(), reltol=1e-8, abstol=1e-10)
sol = integrate(kp, ks)

cell = SemiAdiabaticCell(; Cp=4000.0u"J/K", T_env=293.15u"K", L=0.5u"W/K")
cal  = SemiAdiabaticCalorimeter(; cell=cell, T0=293.15u"K")
sol  = integrate(kp, ks; calorimeter=cal)
t, T_vec = temperature_profile(sol, cal)
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

    # ── Warn for missing Cp° when semi-adiabatic calorimeter requested ─────
    if calorimeter isa SemiAdiabaticCalorimeter
        missing_cp = String[]
        for (sp, cp_fn) in zip(kp.system.species, p.cp_fns)
            isnothing(cp_fn) && push!(missing_cp, string(symbol(sp)))
        end
        if !isempty(missing_cp)
            shown = join(missing_cp[1:min(5, length(missing_cp))], ", ")
            suffix = length(missing_cp) > 5 ? "…" : ""
            @warn "SemiAdiabaticCalorimeter: variable Cp_total requires Cp° data per " *
                "species. Missing for $(length(missing_cp)) species " *
                "($shown$suffix). Their contribution to Cp_total is treated as zero."
        end
    end

    # ── Optionally extend state vector for calorimeter ─────────────────────
    u0 = isnothing(calorimeter) ? u0_base : extend_u0(u0_base, calorimeter)

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
    # H(t) = Σᵢ nᵢ(t)·hᵢ(T) — used by cumulative_heat to compute Q = H(0) - H(t).
    # For SemiAdiabaticCalorimeter the current temperature is u[n_kin+1]; otherwise p_.T.
    is_semi_adiabatic = calorimeter isa SemiAdiabaticCalorimeter
    function _enthalpy_affect!(integrator)
        p_ = integrator.p
        any(!isnothing, p_.h_fns) || return
        for (j, idx) in enumerate(p_.idx_kin_unique)
            p_.n_full[idx] = max(integrator.u[j], p_.ϵ)
        end
        T_curr = is_semi_adiabatic ? integrator.u[n_kin + 1] : p_.T
        H = _total_enthalpy(p_.n_full, p_.h_fns, T_curr)
        push!(p_.saved_H, H)
        push!(p_.saved_t, integrator.t)
        return nothing
    end

    enthalpy_cb = DiscreteCallback(
        Returns(true),
        _enthalpy_affect!;
        initialize = (_, _, _, integrator) -> _enthalpy_affect!(integrator),
        save_positions = (false, false),
    )

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
