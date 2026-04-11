using Pkg
Pkg.activate(@__DIR__)

if !haskey(Pkg.project().dependencies, "OptimaSolver")
    Pkg.develop(path = joinpath(@__DIR__, "..", "..", "OptimaSolver.jl"))
end

using Revise
using ChemistryLab
using Optimization, OptimizationIpopt
using DynamicQuantities

substances_inorg = build_species("data/slop98-inorganic-thermofun.json")
dict = Dict(symbol(s) => s for s in substances_inorg)
species = [dict[s] for s in split("H2O@ Na+ NaOH@ H+ OH-")]

cs = ChemicalSystem(species)

# ── Initial state ──────────────────────────────────────────────────────────────
state = ChemicalState(cs)

set_quantity!(state, "NaOH@", 1.e-3 * u"mol")
set_quantity!(state, "H2O@", 1 * u"kg")
V = volume(state)
set_quantity!(state, "H+",  1e-7u"mol/L" * V.liquid)
set_quantity!(state, "OH-", 1e-7u"mol/L" * V.liquid)
# rescale!(state, 1.0u"kg")

# ── Equilibrate ─────────────────────────────────────────────────────────────────
state_eq = equilibrate(state; variable_space=Val(:linear), verbose=5)

display(state_eq)
