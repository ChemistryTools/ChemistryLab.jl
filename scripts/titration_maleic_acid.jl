using ChemistryLab
using DynamicQuantities
using ProgressMeter

substances_inorg = build_species("data/slop98-inorganic-thermofun.json")
substances_org   = build_species("data/slop98-organic-thermofun.json")

aq_species  = speciation(substances_inorg, split("H2O@ Na+ NaOH@ H+ OH-"); aggregate_state = [AS_AQUEOUS])
ace_species = speciation(substances_org,   split("MalH2@ MalH- Mal-2");             aggregate_state = [AS_AQUEOUS])
species     = unique(s -> symbol(s), vcat(aq_species, ace_species))

dict_all_species = merge(Dict(symbol(s) => s for s in substances_inorg), Dict(symbol(s) => s for s in substances_org))
species = [dict_all_species[s] for s in split("H2O@ Na+ NaOH@ H+ OH- MalH2@ MalH- Mal-2")]

cs = ChemicalSystem(species, ["H2O@", "H+", "Mal-2", "Na+", "Zz"])

using OptimizationIpopt

solver = EquilibriumSolver(
    cs,
    DiluteSolutionModel(),
    IpoptOptimizer(
        # acceptable_tol        = 1e-10,
        # dual_inf_tol          = 1e-10,
        # acceptable_iter       = 1000,
        # constr_viol_tol       = 1e-10,
        mu_strategy = "adaptive",
        # warm_start_init_point = "no",
    );
    variable_space = Val(:linear),
    abstol  = 1e-8,
    reltol  = 1e-8,
    maxiters = 100,
    verbose = 0,
)

V_acid = 100e-3   # volume of acid solution, L
c_acid = 0.1     # maleic acid concentration, mol/L
c_base = 2     # NaOH concentration, mol/L
n_H2A  = V_acid * c_acid   # total moles of H₂A = 2.5 mmol

volumes_NaOH = range(0, 15; length = 201)   # mL
pH_vals = Float64[]

s = ChemicalState(cs)
@showprogress for V_mL in volumes_NaOH
    # @show V_mL
    V_NaOH  = V_mL * 1e-3           # L
    n_NaOH    = c_base * V_NaOH        # mol of NaOH (= mol of Na⁺ added)
    V_total = V_acid + V_NaOH        # total volume, L

    set_quantity!(s, "MalH2@", n_H2A   * u"mol")
    set_quantity!(s, "NaOH@", n_NaOH * u"mol")
    set_quantity!(s, "H2O@",   V_total * u"kg")

    V_liq = volume(s).liquid
    set_quantity!(s, "H+",  1e-7u"mol/L" * V_liq)   # pH-neutral seed
    set_quantity!(s, "OH-", 1e-7u"mol/L" * V_liq)

    s_eq = solve(solver, s)
    push!(pH_vals, pH(s_eq))
end

println("pH at V = 0 mL (pure acid)        : ", round(pH_vals[1],  digits = 2))
println("pH at V = 2.5 mL (½ PE₁, ≈ pKa₁): ", round(pH_vals[6], digits = 2))
println("pH at V = 5 mL  (PE₁)            : ", round(pH_vals[11], digits = 2))
println("pH at V = 7.5 mL (½ PE₂, ≈ pKa₂): ", round(pH_vals[16], digits = 2))
println("pH at V = 10 mL  (PE₂)            : ", round(pH_vals[21], digits = 2))
println("pH at V = 15 mL  (excess NaOH)    : ", round(pH_vals[26], digits = 2))

using Plots #hide

pKa1 = 1.92
pKa2 = 6.27
V_eq1 = n_H2A / c_base * 1e3    # first equivalence point  = 25 mL
V_eq2 = 2 * n_H2A / c_base * 1e3  # second equivalence point = 50 mL

p = plot(
    collect(volumes_NaOH), pH_vals;
    xlabel     = "V(NaOH) (mL)",
    ylabel     = "pH",
    label      = "Titration curve",
    linewidth  = 2,
    marker     = :circle,
    markersize = 3,
    color      = :steelblue,
    title      = "Titration of maleic acid (0.1 M) by NaOH (2 M)",
    ylims      = (0, 14),
    legend     = :topleft,
)
vline!(p, [V_eq1]; linestyle = :dash, color = :red,    label = "PE₁ ($(round(V_eq1,digits = 1)) mL)")
vline!(p, [V_eq2]; linestyle = :dash, color = :blue,   label = "PE₂ ($(round(V_eq2,digits = 1)) mL)")
hline!(p, [pKa1];  linestyle = :dot,  color = :orange, label = "pKₐ₁ = $pKa1")
hline!(p, [pKa2];  linestyle = :dot,  color = :green,  label = "pKₐ₂ = $pKa2")
p
