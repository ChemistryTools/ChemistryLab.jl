"""
    KineticsOrdinaryDiffEqExt

Extension activated when `OrdinaryDiffEq` is loaded. Provides the concrete
`ChemistryLab.integrate(::KineticsProblem, ::KineticsSolver; ...)` implementation
using `OrdinaryDiffEq.ODEProblem` and registers `Rodas5P()` as the default solver.

Follows the Leal et al. (2017) formulation with partial equilibrium and
optional thermal coupling (isothermal or semi-adiabatic calorimetry).

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
    _DEFAULT_KINETICS_SOLVER_FACTORY,
    SemiAdiabaticCalorimeter,
    symbol

# ── Concrete integrate implementation ────────────────────────────────────────

"""
    integrate(kp::KineticsProblem, ks::KineticsSolver) -> ODESolution

Integrate the kinetics ODE using `OrdinaryDiffEq` (Leal et al. 2017 formulation).

The ODE function, initial state, and parameters are built from `kp`.
Calorimetry (isothermal or semi-adiabatic) is integrated directly in the ODE
right-hand-side — no separate `extend_ode!` step.

Default tolerances: `reltol = 1e-8`, `abstol = 1e-10`.

# Examples

```julia
using ChemistryLab, OrdinaryDiffEq

ks  = KineticsSolver(; ode_solver=Rodas5P(), reltol=1e-8, abstol=1e-10)
sol = integrate(kp, ks)
```
"""
function integrate(kp::KineticsProblem, ks::KineticsSolver)
    f! = build_kinetics_ode(kp)
    u0 = build_u0(kp)
    p = build_kinetics_params(kp)

    # Warn for missing Cp° when semi-adiabatic
    if kp.calorimeter isa SemiAdiabaticCalorimeter
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

    defaults = (reltol = 1.0e-8, abstol = 1.0e-10)
    merged = merge(defaults, ks.kwargs)
    solver = isnothing(ks.ode_solver) ? Rodas5P() : ks.ode_solver

    prob = ODEProblem(f!, u0, kp.tspan, p)
    return solve(prob, solver; merged...)
end

# ── __init__: register default solver ────────────────────────────────────────

function __init__()
    return _DEFAULT_KINETICS_SOLVER_FACTORY[] =
        () -> KineticsSolver(; ode_solver = Rodas5P())
end

end  # module KineticsOrdinaryDiffEqExt
