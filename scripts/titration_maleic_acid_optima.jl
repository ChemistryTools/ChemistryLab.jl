using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))   # active l'environnement ChemistryLab.jl

if !haskey(Pkg.project().dependencies, "Optima")
    Pkg.develop(path = joinpath(@__DIR__, "..", "..", "Optima.jl"))
end

using ChemistryLab
using DynamicQuantities
using ProgressMeter
using SciMLBase: solve   # needed: Optima does not re-export SciMLBase.solve
using Optima

substances_inorg = build_species("data/slop98-inorganic-thermofun.json")
substances_org = build_species("data/slop98-organic-thermofun.json")

dict_all_species = merge(
    Dict(symbol(s) => s for s in substances_inorg),
    Dict(symbol(s) => s for s in substances_org),
)
species = [dict_all_species[s] for s in split("H2O@ Na+ NaOH@ H+ OH- MalH2@ MalH- Mal-2")]

cs = ChemicalSystem(species, ["H2O@", "H+", "Mal-2", "Na+", "Zz"])

# ── Solver : Optima with warm-start ────────────────────────────────────────
solver = EquilibriumSolver(
    cs,
    DiluteSolutionModel(),
    OptimaOptimizer(tol = 1.0e-12, warm_start = true);
    variable_space = Val(:linear),
)

V_acid = 100.0e-3  # volume of acid solution, L
c_acid = 0.1     # maleic acid concentration, mol/L
c_base = 2.0     # NaOH concentration, mol/L
n_H2A = V_acid * c_acid  # total moles of H₂A

V_eq1 = n_H2A / c_base * 1.0e3      # first equivalence point,  mL
V_eq2 = 2 * n_H2A / c_base * 1.0e3  # second equivalence point, mL

volumes_NaOH = range(0, 15; length = 201)  # mL
pH_vals = Float64[]

s = ChemicalState(cs)
@showprogress for V_mL in volumes_NaOH
    V_NaOH = V_mL * 1.0e-3           # L
    n_NaOH = c_base * V_NaOH       # mol of NaOH added
    V_total = V_acid + V_NaOH       # total volume, L

    set_quantity!(s, "MalH2@", n_H2A * u"mol")
    set_quantity!(s, "NaOH@", n_NaOH * u"mol")
    set_quantity!(s, "H2O@", V_total * u"kg")

    V_liq = volume(s).liquid
    set_quantity!(s, "H+", 1.0e-7u"mol/L" * V_liq)
    set_quantity!(s, "OH-", 1.0e-7u"mol/L" * V_liq)

    # Near each equivalence point the pH changes steeply; reset to avoid
    # stale warm-start guesses that slow convergence or degrade accuracy.
    # if abs(V_mL - V_eq1) < 0.5 || abs(V_mL - V_eq2) < 0.5
    #     reset_cache!(solver.solver)
    # end

    s_eq = solve(solver, s)
    push!(pH_vals, pH(s_eq))
end

println("pH at V = 0 mL   (pure acid)         : ", round(pH_vals[1], digits = 2))
println("pH at V = 2.5 mL (½ PE₁, ≈ pKa₁)   : ", round(pH_vals[6], digits = 2))
println("pH at V = 5 mL   (PE₁)              : ", round(pH_vals[11], digits = 2))
println("pH at V = 7.5 mL (½ PE₂, ≈ pKa₂)   : ", round(pH_vals[16], digits = 2))
println("pH at V = 10 mL  (PE₂)              : ", round(pH_vals[21], digits = 2))
println("pH at V = 15 mL  (excess NaOH)       : ", round(pH_vals[26], digits = 2))

using Plots

pKa1 = 1.92
pKa2 = 6.27

plt = plot(
    collect(volumes_NaOH), pH_vals;
    xlabel = "V(NaOH) (mL)",
    ylabel = "pH",
    label = "Titration curve (Optima)",
    linewidth = 2,
    marker = :circle,
    markersize = 3,
    color = :steelblue,
    title = "Titration of maleic acid (0.1 M) by NaOH (2 M)",
    ylims = (0, 14),
    legend = :topleft,
)
vline!(plt, [V_eq1]; linestyle = :dash, color = :red, label = "PE₁ ($(round(V_eq1, digits = 1)) mL)")
vline!(plt, [V_eq2]; linestyle = :dash, color = :blue, label = "PE₂ ($(round(V_eq2, digits = 1)) mL)")
hline!(plt, [pKa1]; linestyle = :dot, color = :orange, label = "pKₐ₁ = $pKa1")
hline!(plt, [pKa2]; linestyle = :dot, color = :green, label = "pKₐ₂ = $pKa2")
plt
