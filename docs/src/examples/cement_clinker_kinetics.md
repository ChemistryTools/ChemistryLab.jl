# Cement clinker hydration kinetics

This example demonstrates the full kinetics workflow: from database loading to
ODE integration and calorimetric post-processing. It models the hydration of an
OPC (CEM I 52.5 R) clinker using the Parrot--Killoh [ParrotKilloh1984](@cite)
rate model with Arrhenius temperature correction
[SchindlerFolliard2005](@cite), coupled to a semi-adiabatic calorimeter
[Lavergne2018](@cite).

## 1. Chemical system

We select the four clinker phases as kinetic species, plus the main hydration
products and water. `CEMDATA_PRIMARIES` provides the independent aqueous
components.

```@example clinker
using ChemistryLab
using OrdinaryDiffEq
using DynamicQuantities
using OrderedCollections

DATA_FILE = joinpath(pkgdir(ChemistryLab), "data", "cemdata18-thermofun.json")

substances = build_species(DATA_FILE)

input_species = split(
    "C3S C2S C3A C4AF " *
        "Portlandite Jennite ettringite monosulphate12 C3AH6 C3FH6 " *
        "H2O@",
)

species = speciation(substances, input_species; aggregate_state = [AS_AQUEOUS])
cs = ChemicalSystem(species, CEMDATA_PRIMARIES)

println("Chemical system: $(length(cs.species)) species")
```

## 2. Initial state

Typical CEM I 52.5 R composition from [Lavergne2018](@cite). We work with 1 kg
of cement at water/cement ratio ``w/c = 0.4``.

```@example clinker
const WC = 0.4
const COMPOSITION = (C3S = 0.619, C2S = 0.165, C3A = 0.08, C4AF = 0.087)

state0 = ChemicalState(cs)
for (name, frac) in pairs(COMPOSITION)
    set_quantity!(state0, string(name), frac * u"kg")
end
set_quantity!(state0, "H2O@", WC * u"kg")
nothing # hide
```

## 3. Parrot--Killoh kinetic models

Maximum degree of hydration from the [Powers1948](@cite) law:
``\alpha_{\max} \leq w/c \,/\, 0.42``.

```@example clinker
const α_max = min(1.0, WC / 0.42)

pk_C3S  = parrot_killoh(PK_PARAMS_C3S,  "C3S";  α_max)
pk_C2S  = parrot_killoh(PK_PARAMS_C2S,  "C2S";  α_max)
pk_C3A  = parrot_killoh(PK_PARAMS_C3A,  "C3A";  α_max)
pk_C4AF = parrot_killoh(PK_PARAMS_C4AF, "C4AF"; α_max)
nothing # hide
```

## 4. Balanced hydration reactions

Reactions must be **mass-balanced** so that ``\Delta_r H^0`` can be computed
automatically from species ``\Delta_a H^0``. The Jennite formula unit in
CEMDATA18 is ``(\text{SiO}_2)(\text{CaO})_{5/3}(\text{H}_2\text{O})_{21/10}``.

```@example clinker
sp(name) = cs[name]

rxn_C3S = Reaction(
    OrderedDict(sp("C3S") => 1.0, sp("H2O@") => 103/30),
    OrderedDict(sp("Jennite") => 1.0, sp("Portlandite") => 4/3);
    symbol = "C₃S hydration",
)
rxn_C3S[:rate] = pk_C3S

rxn_C2S = Reaction(
    OrderedDict(sp("C2S") => 1.0, sp("H2O@") => 73/30),
    OrderedDict(sp("Jennite") => 1.0, sp("Portlandite") => 1/3);
    symbol = "C₂S hydration",
)
rxn_C2S[:rate] = pk_C2S

rxn_C3A = Reaction(
    OrderedDict(sp("C3A") => 1.0, sp("H2O@") => 6.0),
    OrderedDict(sp("C3AH6") => 1.0);
    symbol = "C₃A hydration",
)
rxn_C3A[:rate] = pk_C3A

rxn_C4AF = Reaction(
    OrderedDict(sp("C4AF") => 1.0, sp("Portlandite") => 2.0, sp("H2O@") => 10.0),
    OrderedDict(sp("C3AH6") => 1.0, sp("C3FH6") => 1.0);
    symbol = "C₄AF hydration",
)
rxn_C4AF[:rate] = pk_C4AF

kinetic_reactions = [rxn_C3S, rxn_C2S, rxn_C3A, rxn_C4AF]
nothing # hide
```

## 5. Kinetics problem and calorimeter

The semi-adiabatic calorimeter [Lavergne2018](@cite) solves
```math
\frac{dT}{dt} = \frac{\dot{q}(t) - \varphi(\Delta T)}{C_p^{\text{total}}}
```
where ``\varphi(\Delta T) = a \, \Delta T + b \, \Delta T^2`` models quadratic
heat losses and ``C_p^{\text{total}} = C_p + \sum_i n_i C_{p,i}^0(T)``.

```@example clinker
cal = SemiAdiabaticCalorimeter(;
    Cp        = (1.0 * 800.0 + WC * 4186.0 + 1.0 * 900.0) * u"J/K",
    T_env     = 293.15u"K",
    heat_loss = ΔT -> 0.3 * ΔT + 0.003 * ΔT^2,
    T0        = 293.15u"K",
)

kp = KineticsProblem(
    cs, kinetic_reactions, state0, (0.0, 7.0 * 86400.0);
    calorimeter = cal,
    equilibrium_solver = nothing,
)
nothing # hide
```

## 6. Integration

```@example clinker
ks  = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1e-6, abstol = 1e-9)
sol = integrate(kp, ks)
println("$(length(sol.t)) accepted steps over 7 days")
```

## 7. Post-processing

```@example clinker
using Printf

t_h = sol.t ./ 3600.0

t_T, T_K_vec = temperature_profile(sol, cal)
t_Q, Q_J_vec = cumulative_heat(sol, cal)
T_°C = T_K_vec .- 273.15
Q_kJ = Q_J_vec ./ 1000.0

n0_kin = [sol.prob.p.n_initial_full[i] for i in kp.idx_kinetic]
n_kin  = [[u[i] for u in sol.u] for i in eachindex(n0_kin)]

function phase_alpha(cs, kp, sol, n0_kin, n_kin, name)
    sp_idx = findfirst(sp -> ChemistryLab.symbol(sp) == name, cs.species)
    pos = findfirst(==(sp_idx), kp.idx_kinetic)
    isnothing(pos) && return fill(NaN, length(sol.t))
    return 1.0 .- n_kin[pos] ./ n0_kin[pos]
end

α_C3S  = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C3S")
α_C2S  = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C2S")
α_C3A  = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C3A")
α_C4AF = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C4AF")

w = COMPOSITION
α_mean = (w.C3S .* α_C3S .+ w.C2S .* α_C2S .+
          w.C3A .* α_C3A .+ w.C4AF .* α_C4AF) ./
         (w.C3S + w.C2S + w.C3A + w.C4AF)

@printf "ΔT max   = %.2f °C\n" maximum(T_°C) - 20.0
@printf "α(C₃S)   = %.4f\n"   α_C3S[end]
@printf "α(C₂S)   = %.4f\n"   α_C2S[end]
@printf "α(C₃A)   = %.4f\n"   α_C3A[end]
@printf "α(C₄AF)  = %.4f\n"   α_C4AF[end]
@printf "ᾱ mean   = %.4f\n"   α_mean[end]
@printf "Q total  = %.1f kJ/kg\n" Q_kJ[end]
```

## 8. Plots

```@example clinker
using Plots
gr()

p1 = plot(t_T ./ 3600, T_°C;
    xlabel="Time [h]", ylabel="T [°C]",
    title="Temperature", label="T(t)", lw=2, color=:red)
hline!(p1, [20.0]; ls=:dash, color=:gray, label="T₀")

p2 = plot(t_h, [α_C3S α_C2S α_C3A α_C4AF α_mean];
    xlabel="Time [h]", ylabel="α",
    title="Degree of hydration", lw=2,
    label=["C₃S" "C₂S" "C₃A" "C₄AF" "ᾱ"],
    ls=[:solid :dash :dot :dashdot :solid])
hline!(p2, [α_max]; ls=:dash, color=:black, label="α_max")

p3 = plot(t_Q ./ 3600, Q_kJ;
    xlabel="Time [h]", ylabel="Q [kJ/kg]",
    title="Cumulative heat", label="Q(t)", lw=2, color=:purple)

plot(p1, p2, p3; layout=(1,3), size=(1400, 420),
    plot_title="CEM I w/c=$WC — Parrot–Killoh + semi-adiabatic calorimeter")
```
