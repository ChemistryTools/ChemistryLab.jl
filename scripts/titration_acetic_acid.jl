using ChemistryLab
using DynamicQuantities

substances_inorg = build_species("data/slop98-inorganic-thermofun.json")
substances_org   = build_species("data/slop98-organic-thermofun.json")

dict_all_species = merge(Dict(symbol(s) => s for s in substances_inorg), Dict(symbol(s) => s for s in substances_org))
species = [dict_all_species[s] for s in split("H2O@ Na+ NaOH@ H+ OH- AceH@ Ace-")]

cs = ChemicalSystem(species, ["H2O@", "H+", "Ace-", "Na+"])

R  = ustrip(Constants.R)
T  = 298.15   # K
RT = R * T

sp = Dict(symbol(s) => s for s in cs.species)

Ka  = exp(-(sp["Ace-"].ΔₐG⁰(T = T) + sp["H+"].ΔₐG⁰(T = T) - sp["AceH@"].ΔₐG⁰(T = T)) / RT)
pKa = -log10(Ka)
println("Ka  = ", Ka)
println("pKa = ", round(pKa, digits = 2))

using OptimizationIpopt #hide

solver = EquilibriumSolver(
    cs,
    DiluteSolutionModel(),
    IpoptOptimizer(
        acceptable_tol        = 1e-10,
        dual_inf_tol          = 1e-10,
        acceptable_iter       = 1000,
        constr_viol_tol       = 1e-10,
        warm_start_init_point = "no",
    );
    variable_space = Val(:linear),
    abstol  = 1e-10,
    reltol  = 1e-10,
    verbose = 0,
)

ca   = 0.1     # acetic acid concentration, mol/L
Va   = 0.1     # volume of acid solution, L
cb   = 0.1     # NaOH concentration, mol/L
nAH  = ca * Va # total moles of CH₃COOH = 10 mmol

Vbeq = nAH / cb          # equivalence volume, L
V_eq = Vbeq * 1e3        # equivalence volume, mL

volumes_NaOH = range(0, 2 * V_eq; length = 100)   # mL
pH_vals = Float64[]

s = ChemicalState(cs)
for V_mL in volumes_NaOH
    # @show V_mL
    Vb      = V_mL * 1e-3     # L
    n_NaOH  = cb * Vb         # mol of NaOH added
    V_total = Va + Vb          # total volume, L

    set_quantity!(s, "AceH@", nAH    * u"mol")
    set_quantity!(s, "NaOH@", n_NaOH * u"mol")
    set_quantity!(s, "H2O@",  V_total * u"kg")

    V_liq = volume(s).liquid
    set_quantity!(s, "H+",  1e-7u"mol/L" * V_liq)   # pH-neutral seed
    set_quantity!(s, "OH-", 1e-7u"mol/L" * V_liq)

    s_eq = solve(solver, s)
    push!(pH_vals, pH(s_eq))
end

idx_heq = argmin(abs.(collect(volumes_NaOH) .- V_eq / 2))
idx_eq  = argmin(abs.(collect(volumes_NaOH) .- V_eq))

println("pH at V = 0 mL    (pure acid)        : ", round(pH_vals[begin],   digits = 2))
println("pH at V = 50 mL   (½ PE, ≈ pKₐ)     : ", round(pH_vals[idx_heq], digits = 2))
println("pH at V = 100 mL  (equivalence point): ", round(pH_vals[idx_eq],   digits = 2))
println("pH at V = 300 mL  (excess NaOH)      : ", round(pH_vals[end],      digits = 2))

using Plots #hide

function ana_pH(Vb_mL)
    Vb = Vb_mL * 1e-3
    Vb ≈ 0 && return (-log10(ca) + pKa) / 2   # pure acid approximation
    return Vb < Vbeq ?
        pKa + log10(cb * Vb / (nAH - cb * Vb)) :
        14  + log10((cb * Vb - nAH) / (Va + Vb))
end

pH0  = (-log10(ca) + pKa) / 2
pHeq = (pKa + 14 + log10(ca * Va / (Va + Vbeq))) / 2

p = plot(
    collect(volumes_NaOH), pH_vals;
    xlabel     = "V(NaOH) (mL)",
    ylabel     = "pH",
    label      = "Numerical (ChemistryLab)",
    linewidth  = 0,
    marker     = :circle,
    markersize = 3,
    color      = :steelblue,
    title      = "Titration of acetic acid (0.1 M) by NaOH (0.1 M)",
    ylims      = (0, 14),
    legend     = :bottomright,
)
lVb = range(0.01 * V_eq, 3 * V_eq; length = 10_000)
plot!(p, lVb, ana_pH.(lVb);
    label     = "Analytical (Henderson–Hasselbalch)",
    linewidth = 2,
    color     = :orange,
)
scatter!(p, [0, V_eq], [pH0, pHeq];
    label      = "Characteristic points",
    color      = :red,
    markersize = 6,
)
vline!(p, [V_eq]; linestyle = :dash, color = :red,  label = "PE ($(round(V_eq,digits = 1)) mL)")
hline!(p, [pKa];  linestyle = :dot,  color = :grey, label = "pKₐ = $(round(pKa, digits = 2)) at V = 50 mL")
p
