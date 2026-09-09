# [Effect of Water/Cement Ratio on Cement Hydration](@id sec-wc-ratio)

The **water-to-cement ratio** (w/c) is the single most important mix-design parameter
of concrete. It controls workability, compressive strength, and durability simultaneously.
From a thermodynamic perspective, it determines how much water is available to hydrate the
clinker phases, which in turn governs the nature and amount of hydration products, the
porosity of the paste, and the pH of the pore solution.

This example scans w/c from 0.30 to 0.60 and tracks, at **full thermodynamic
equilibrium**, the pH of the pore solution, the hydrate assemblage, and the
porosity of the hardened paste in the sealed-curing convention.

!!! warning "Equilibrium answers a narrower question than mix design does"
    Read the scan below for what it is: the assemblage a paste would reach if
    every reaction ran to completion. It is **not** the state of a real paste at
    an age, and two of the effects usually attributed to w/c are absent from it
    by construction — over the range scanned no clinker survives, and there is no
    optimum. Those are **kinetic** limitations, carried by
    [`powers_alpha_max`](@ref) and the rate laws of the
    [kinetics tutorial](@ref sec-kinetics), not by the Gibbs minimum. A
    water-limited regime does exist in the Gibbs minimum, but it starts far below
    Powers' 0.42 — measured on this species list, between w/c = 0.28 and 0.30 —
    and the two limits are not the same statement, which the analysis below
    takes apart.
    This page is the reference the kinetic calculation converges toward; the
    [coupled hydration example](@ref sec-coupled-hydration) is the one to read
    for an age.

---

## System setup

The same clinker composition and species set as the "simplified clinker dissolution" example are used here.

```@example wc_setup
using ChemistryLab
using DynamicQuantities

substances = build_species(datapath("cemdata18-thermofun.json"))

input_species = split("C3S C2S C3A C4AF Gp Anh Portlandite Jennite H2O@ ettringite monosulphate12 C3AH6 C3FH6 C4FH13")
species = speciation(substances, input_species; aggregate_state = [AS_AQUEOUS])

cs = ChemicalSystem(species, CEMDATA_PRIMARIES)
```

The clinker composition (mass fractions of the anhydrous cement phases) is fixed throughout the scan:

| Phase | Symbol | Mass fraction |
|:------|:-------|:--------------|
| Alite | `C3S`  | 67.8 % |
| Belite | `C2S` | 16.6 % |
| Aluminate | `C3A` | 4.0 % |
| Ferrite | `C4AF` | 7.2 % |
| Gypsum | `Gp`  | 2.8 % |

```@example wc_setup
compo = ["C3S" => 0.678, "C2S" => 0.166, "C3A" => 0.040, "C4AF" => 0.072, "Gp" => 0.028]
c     = sum(last.(compo))   # cement mass fraction (= 0.984 here)
```

---

## The solver

[`equilibrate_certified`](@ref) is used throughout, which needs nothing but
`OptimaSolver` loaded:

```@example wc_setup
using OptimaSolver
```

That route returns `(state, certificate)`, and the certificate is what makes the
numbers below quotable: for a convex problem the KKT conditions are *sufficient*,
so `cert.optimal == true` is a **proof** that the composition is the Gibbs
minimum and not merely the point an iteration stopped at. Every value on this
page is checked that way, and the check is not idle — on a cement the difference
between "the solver returned" and "the answer is proved" has been measured at
13 % of the total volume and four units of pH.

!!! note "Solving through Ipopt instead"
    An `EquilibriumSolver` around any nonlinear back end still works, and reaches
    the same assemblage here:

    ```julia
    using Optimization, OptimizationIpopt

    opt = IpoptOptimizer(
        acceptable_tol        = 1e-10,
        dual_inf_tol          = 1e-10,
        acceptable_iter       = 100,
        constr_viol_tol       = 1e-10,
        warm_start_init_point = "no",
    )
    solver = EquilibriumSolver(
        cs, DiluteSolutionModel(), opt;
        variable_space = Val(:linear), abstol = 1e-8, reltol = 1e-8,
    )
    eq = solve(solver, deepcopy(fresh))     # no certificate
    ```

    It is left out of the executed page for two reasons, neither of them about
    the quality of Ipopt as an optimizer. It is a **bare interior point**: it
    returns an iterate and no statement about it, so a caller has to audit the
    answer with [`optimality_certificate`](@ref) anyway. And it is an extra
    binary dependency for a calculation the package can already prove. The
    visible difference on this page is small but telling: the interior point
    leaves the absent phases at its lower bound, around `1e-8` mol, while the
    certified route puts them at exactly zero.

---

## Scanning the w/c ratio

For each value of w/c the fresh state is rebuilt from scratch, and the total mass is
normalized to 1 kg of paste (cement + water) so that all amounts are comparable
across the scan. The fresh state is kept: it is the volume reference the porosity
is referred to.

```@example wc_setup
sp_idx   = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))
wc_range = range(0.30, 0.60; length = 13)

function fresh_paste(wc)
    w, mtot = wc * c, c + wc * c
    st = ChemicalState(cs)
    for (sym, mfrac) in compo
        set_quantity!(st, sym, mfrac / mtot * u"kg")
    end
    set_quantity!(st, "H2O@", w / mtot * u"kg")
    V = volume(st)
    set_quantity!(st, "H+",  1e-7u"mol/L" * V.liquid)   # charge seed, pH-neutral
    set_quantity!(st, "OH-", 1e-7u"mol/L" * V.liquid)
    return st
end

pH_vals    = Float64[]
ϕ_liquid   = Float64[]
ϕ_void     = Float64[]
ϕ_total    = Float64[]
n_portl    = Float64[]
n_mono     = Float64[]
n_ett      = Float64[]
n_jennite  = Float64[]
n_clinker  = Float64[]

certified = Bool[]

for wc in wc_range
    fresh    = fresh_paste(wc)
    eq, cert = equilibrate_certified(deepcopy(fresh))   # the solve may mutate its argument
    push!(certified, cert.optimal)
    ϕ = porosity(eq, fresh)                  # sealed-curing convention
    amount(sym) = ustrip(eq.n[sp_idx[sym]])
    push!(pH_vals,   pH(eq))
    push!(ϕ_liquid,  ϕ.liquid)
    push!(ϕ_void,    ϕ.void)
    push!(ϕ_total,   ϕ.total)
    push!(n_portl,   amount("Portlandite"))
    push!(n_mono,    amount("monosulphate12"))
    push!(n_ett,     amount("ettringite"))
    push!(n_jennite, amount("Jennite"))
    push!(n_clinker, sum(amount(s) for s in ("C3S", "C2S", "C3A", "C4AF")))
end
nothing # hide
```

!!! note "Why `porosity(eq, fresh)` and not `porosity(eq)`"
    The one-argument [`porosity`](@ref) divides the pore volume by the **current**
    total volume, which is right for a fixed-volume aqueous system and wrong for a
    setting binder: hydration products occupy less space than the reactants they
    consume, so the paste's own volume shrinks under the ratio. The two-argument
    form refers everything to the fresh volume and splits the result into the
    water-filled porosity and the empty porosity that the Le Chatelier
    contraction creates. On this mix at w/c = 0.50 the one-argument form returns
    0.285 against a total of 0.338, a gap of 5.3 points of porosity, the total
    volume having shrunk by 7.4 %.

---

## Results

### pH and porosity

```@example wc_setup
using Plots

p1 = plot(collect(wc_range), pH_vals;
    xlabel = "w/c ratio", ylabel = "Pore solution pH", label = "pH",
    linewidth = 2, marker = :circle, markersize = 4, color = :steelblue,
    title = "Pore solution pH", ylims = (11.5, 13.5), legend = :bottomright)

p2 = plot(collect(wc_range), ϕ_total .* 100;
    xlabel = "w/c ratio", ylabel = "Porosity (%)", label = "total",
    linewidth = 2, marker = :circle, markersize = 4, color = :firebrick,
    title = "Porosity referred to the fresh volume",
    ylims = (0, 50), legend = :topleft)
plot!(p2, collect(wc_range), ϕ_liquid .* 100;
    label = "water-filled", linewidth = 2, marker = :square, markersize = 3,
    color = :steelblue)
plot!(p2, collect(wc_range), ϕ_void .* 100;
    label = "empty (Le Chatelier)", linewidth = 2, marker = :diamond,
    markersize = 3, color = :seagreen)

plot(p1, p2; layout = (1, 2), left_margin = 8Plots.mm,
     bottom_margin = 8Plots.mm, size = (950, 410))
```

The pH is **flat to three decimals** over the whole range. That is the signature
of a buffered solution: portlandite is present at every w/c, so the calcium and
hydroxide activities are pinned by its saturation, and in a dilute-solution model
those activities do not know how large the pore volume is. Diluting a saturated
solution with more of its own solvent does not change its pH — it dissolves more
portlandite. The pH would start to move only once portlandite is exhausted, which
this composition never does.

### Phase assemblage

```@example wc_setup
p3 = plot(collect(wc_range), n_portl;
    xlabel = "w/c ratio", ylabel = "Amount (mol / kg of paste)",
    label = "Portlandite  Ca(OH)₂", linewidth = 2, marker = :circle,
    markersize = 4, color = :steelblue, title = "Hydrate assemblage at equilibrium",
    legend = :right)
plot!(p3, collect(wc_range), n_jennite;
    label = "Jennite (C-S-H)", linewidth = 2, marker = :square,
    markersize = 3, color = :firebrick)
plot!(p3, collect(wc_range), n_mono;
    label = "monosulphate12 (AFm)", linewidth = 2, marker = :diamond,
    markersize = 3, color = :seagreen)
plot!(p3, collect(wc_range), n_ett;
    label = "ettringite (AFt)", linewidth = 2, marker = :utriangle,
    markersize = 3, color = :darkorange)
plot(p3; left_margin = 8Plots.mm, bottom_margin = 8Plots.mm, size = (700, 420))
```

```@example wc_setup
using Printf
@printf "every point certified                    : %s\n" all(certified)
@printf "ettringite, largest value over the scan  : %.3e mol\n" maximum(n_ett)
@printf "clinker left, largest value over the scan: %.3e mol\n" maximum(n_clinker)
```

---

## Analysis

| Quantity | What the scan gives | Why |
|:--|:--|:--|
| **pH** | 12.39, constant to three decimals | buffered by portlandite saturation; independent of pore volume in a dilute model |
| **Portlandite, C-S-H, AFm** | decrease with w/c, in proportion to the cement fraction | amounts are per kg of *paste*; adding water dilutes the binder, it does not change what a gram of cement produces |
| **Ettringite** | exactly zero at every w/c | see below |
| **Clinker left** | exactly zero at every w/c of *this* scan | see below, and it is not "at any w/c" |
| **Total porosity** | 12.2 % to 41.0 %, monotone, no optimum | the excess water has nowhere to go but the pore space |
| **Empty porosity** | 9.8 % down to 6.6 % | the Le Chatelier contraction is roughly fixed per gram of cement, so it is a smaller fraction of a larger reference volume |

Two of those deserve stating plainly, because the intuition of mix design points
the other way.

!!! warning "Ettringite does not form here, and that is correct"
    AFt needs about three sulfates per aluminate; this clinker has 2.8 % gypsum
    against 4.0 % C₃A plus the ferrite, so at **equilibrium** the sulfate is all
    taken up by the AFm phase `monosulphate12`, and `ettringite` comes back at
    exactly zero — absence, not a small amount. Ettringite is the phase that forms **early**, while
    sulfate is still locally abundant, and then converts to AFm as it runs out.
    It is a kinetic intermediate for this mix, so an equilibrium scan cannot show
    it. Raise the gypsum content and AFt becomes stable — that is the
    sulfate-balance calculation, and it is worth doing before reading anything
    into an AFt amount.

!!! warning "No clinker survives *over this range*, and Powers' 0.42 is not the reason"
    The four anhydrous phases come back at exactly zero at every w/c of the scan,
    0.30 included, and the certificate proves it. That looks like it contradicts
    [Powers1948](@cite), whose ``\alpha_{\max} = w/c \,/\, 0.42`` leaves
    unreacted clinker in any paste below w/c = 0.42 — 0.36 with curing water. It
    does not, and the reason is worth stating because the two numbers measure
    different things.

    Powers' 0.42 g of water per gram of cement is **not** a stoichiometric
    demand. It is about 0.23 g of *non-evaporable* water, which is the water
    written into the hydrate formulae, plus about 0.19 g of **gel water** held in
    the C-S-H gel pores. Only the first is a mass balance that a Gibbs
    minimization must respect. The second is water that is physically there and
    chemically unavailable: in a sealed paste it is immobilized in pores too fine
    to feed further reaction, and hydration stops by self-desiccation with water
    still in the specimen. That is a statement about **transport and access**,
    which no equilibrium calculation contains, and it is what
    [`powers_alpha_max`](@ref) carries into the kinetic rate laws.

    So between roughly 0.23 and 0.42 the two disagree **and both are right**: the
    water suffices to write the hydrates, and a real sealed paste still cannot
    reach them. Below the stoichiometric demand they agree, because there the
    limit is mass balance and the Gibbs minimum obeys it like anything else.

    Where that crossover falls is a property of the hydrate assemblage, not a
    constant, so it is measured rather than quoted:

```@example wc_setup
# Below the scanned range the water runs out, and the minimum leaves clinker.
low = [0.15, 0.20, 0.25, 0.28, 0.30]
println(" w/c   clinker left (% of the cement)   certified")
for wc in low
    fr = fresh_paste(wc)
    m0 = sum(
        ustrip(us"kg", fr.n[sp_idx[s]] * cs.species[sp_idx[s]][:M])
            for s in ("C3S", "C2S", "C3A", "C4AF")
    )
    # The warnings are the point of the last column, not of the transcript.
    eq, cert = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        equilibrate_certified(deepcopy(fr))
    end
    ml = sum(
        ustrip(us"kg", eq.n[sp_idx[s]] * cs.species[sp_idx[s]][:M])
            for s in ("C3S", "C2S", "C3A", "C4AF")
    )
    @printf("%5.2f   %26.1f   %9s\n", wc, 100 * ml / m0, cert.optimal)
end
```

!!! warning "Read the last column"
    The water-limited points are **not certified**. Their amounts are the best
    answers the multi-start route found, not proved minima, and the reason is the
    one the [equilibrium tutorial](@ref sec-equilibrium) gives for cold starts in
    general: a water-starved paste is the badly conditioned end of this problem.
    What does not depend on the solver is the **existence** of the regime — 5 g
    of water cannot hydrate 100 g of cement whatever the algorithm, since the
    hydrates would need some 23 g. Read the table for that, and for the order of
    magnitude; do not quote its numbers as equilibrium values.

    The previous version of this page said no clinker survives *at any* w/c and
    explained it by the minimum "always forming a less hydrous assemblage". The
    first half is true only over the range scanned, and the second is false: the
    least hydrous assemblage available still binds water, and when there is not
    enough, alite stays.

    What remains true is the rest of the original claim, and it matters for mix
    design: over 0.30–0.60 this scan shows **no optimum w/c and no inflection**.
    The porosity rises monotonically and the minimum-porosity mix design does not
    appear, because it is set by the degree of hydration a paste actually
    reaches, not by the assemblage it would reach given time.

### Assumptions behind these numbers

  - **Complete reaction.** Full equilibrium, no time, no kinetic barrier.
  - **Sealed curing.** No water exchanged with the outside, so the empty
    porosity from the chemical shrinkage stays empty. An immersed specimen would
    draw water in and `ϕ.void` would fill.
  - **The fresh volume is the reference**, and it is held fixed — the specimen
    keeps its cast dimensions, the contraction showing up as internal void rather
    than as shrinkage of the outside. This is the usual convention for a set
    paste; it is wrong before setting, when the material still contracts
    externally.
  - **Ideal molar volumes.** Phase volumes are the sum of ``n_i V_i^0``, with no
    mixing term.
  - **Dilute solution model**, so no ionic-strength correction. The pore solution
    of a cement paste is around 0.1–0.3 mol/kg, where activity coefficients
    depart from unity by tens of percent — use [`HKFActivityModel`](@ref) or
    [`DaviesActivityModel`](@ref) if the ion concentrations themselves matter.
    The pH being buffered, it is the quantity least affected by this choice.
  - **The species list is closed.** Only the 14 species selected above may form;
    siliceous hydrogarnet, hydrotalcite and the alkali sulfates are absent, as
    are the alkalis themselves, which in a real paste raise the pore-solution pH
    to 13 or above.

!!! note "Extending the scan"
    To study **supplementary cementitious materials**, substitute part of the
    clinker and add the corresponding species from `cemdata18` (e.g. `C2ASH8`).
    For **carbonation**, add `CO2@`, `HCO3-`, `CO3-2` and the carbonate phases —
    see the [cement carbonation example](@ref sec-cement-carbonation).
