# [Self-desiccation: where Powers' 0.42 comes from](@id sec-self-desiccation)

A sealed cement paste stops hydrating before it runs out of cement. Powers'
rule of thumb caps the degree of hydration at

```math
\alpha_{\max} = \frac{w/c}{0.42},
```

and [`powers_alpha_max`](@ref) supplies it to the rate laws. The coefficient is
empirical, and this page takes it apart: it writes the water budget that produces
it, fills each term from an **independent** source — the chemistry from a Gibbs
energy minimization, the arrest saturation from a published desorption isotherm —
and reports what the two together imply.

Nothing below is fitted to 0.42. The
[companion script](https://github.com/MicroPoroChemoMechanics/ChemistryLab.jl/blob/main/scripts/self_desiccation_powers.jl)
runs the whole calculation outside the documentation.

!!! warning "Read this first: the Kelvin term does not arrest hydration"
    It is tempting to expect the arrest from thermodynamics — water in a fine
    pore is held at a reduced activity, so hydrates that consume water become
    less stable. **The magnitude is wrong by two orders of magnitude.** At
    ``a_w = 0.80`` the shift is ``RT \ln a_w = -553`` J per mole of water, and
    alite going to a C-S-H plus portlandite consumes about 3.3 mol of water per
    mole, so the whole capillary contribution to the affinity is under 2 kJ/mol
    against a hydration Gibbs energy of order ``-100`` kJ/mol. Nulling it would
    need ``a_w \approx 5\times10^{-6}``, a Kelvin radius smaller than a water
    molecule. Measured in §7 on this page's own paste, with a certificate on
    every answer: imposing a water activity from saturation down to 0.80 leaves
    the equilibrium assemblage **unchanged to six digits**.

    A real paste stops at 75–80 % relative humidity because transport and
    nucleation stop. So the humidity belongs in the **rate law**, through
    [`humidity_factor`](@ref) and [`PoreHumidity`](@ref), and
    [`CapillaryWater`](@ref) is what makes the water activity of the equilibrium
    state mean the pore water rather than the mole fraction. This page uses the
    thermodynamics for the **water budget** and the humidity only as the arrest
    criterion.

---

## 1. The water budget

Take one gram of cement and ``w/c`` grams of mixing water, sealed: no water
enters or leaves. At a degree of hydration ``\alpha``, three things have happened
to that water.

Some of it is **bound into the hydrate formulae** — the ``\ce{H2O}`` of
``\ce{Ca(OH)2}``, of the C-S-H, of ettringite. Write that ``b\,\alpha`` grams per
gram of cement, with ``b`` the demand at full hydration.

What is left is **pore solution**, and its volume follows:

```math
V_\ell = \frac{w/c - b\,\alpha}{\rho_w}.
```

And the products occupy **less** space than the reactants they consumed — the
Le Chatelier contraction. In a sealed specimen that deficit cannot be filled from
outside, so it becomes empty porosity, ``s\,\alpha`` cubic centimeters per gram
of cement, with ``s`` the chemical shrinkage at full hydration.

The pore space is therefore liquid plus void, and its **degree of saturation** is

```math
S(\alpha) = \frac{V_\ell}{V_\ell + V_{\text{void}}}
          = \frac{w/c - b\,\alpha}{w/c - b\,\alpha + s\,\alpha}.
```

Now suppose hydration stops when the saturation falls to some value ``S^\ast``.
Setting ``S(\alpha_{\max}) = S^\ast`` and solving:

```math
(w/c)(1 - S^\ast) = \alpha_{\max}\left[b(1 - S^\ast) + s\,S^\ast\right]
```

```math
\boxed{\;\alpha_{\max} = \frac{w/c}{k}, \qquad
       k = b + s\,\frac{S^\ast}{1 - S^\ast}.\;}
```

Two things follow, and the first is a warning about what this page can prove.

!!! note "Powers' *form* is structural; only the *coefficient* is predicted"
    ``\alpha_{\max} \propto w/c`` came out of the algebra without any physics
    beyond "the budget is linear in ``\alpha``" and "the arrest happens at a
    fixed saturation". Any monotone retention curve and any fixed humidity
    threshold reproduce that proportionality exactly. Finding a straight line is
    therefore **no evidence at all**.

    What is genuinely predicted is ``k``, and ``k`` needs three numbers: ``b``
    and ``s`` from the chemistry, ``S^\ast`` from a measurement. Those are what
    the rest of this page supplies.

---

## 2. The chemistry: ``b`` and ``s``

The cement is the one whose isotherm this page uses further down — taking the
isotherm from one paper and the clinker from another would compare two materials.
[BaroghelBouny1999](@cite) give its mineral composition in their Table 2.

```@example sd
using ChemistryLab
using DynamicQuantities
using OptimaSolver
using Printf

# Baroghel-Bouny et al. (1999), Table 2. The balance is free lime and alkalis,
# which this species list does not carry.
compo = ["C3S" => 0.5728, "C2S" => 0.2398, "C3A" => 0.0303,
         "C4AF" => 0.0759, "Gp" => 0.0439, "Cal" => 0.0184]
cmass = sum(last.(compo))
wc    = 0.34                     # their mix CO
M_H2O = 0.0180153                # kg/mol

phases = split("C3S C2S C3A C4AF Gp Anh Cal Portlandite Jennite H2O@ " *
               "ettringite monosulphate12 C3AH6 C3FH6 C4FH13 monocarbonate")
substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
cs = ChemicalSystem(
    speciation(substances, phases; aggregate_state = [AS_AQUEOUS]), CEMDATA_PRIMARIES
)
iw = only(cs.idx_solvent)
ic = [findfirst(s -> symbol(s) == sym, cs.species) for (sym, _) in compo]
nothing # hide
```

An equilibrium calculation has no notion of a reaction that has not happened, so
the degree of hydration is **imposed**: a fraction ``\alpha`` of the cement is
made available to react with *all* of the mixing water, and the unreacted
remainder is added back afterwards for the volume balance. This is the
construction of [LothenbachWinnefeld2006](@cite), and it is the same one the
[w/c example](@ref sec-wc-ratio) uses below its stoichiometric water demand.

```@example sd
function paste(α)
    mtot = cmass + wc * cmass
    st = ChemicalState(cs)
    for (sym, mfrac) in compo
        set_quantity!(st, sym, α * mfrac / mtot * u"kg")
    end
    set_quantity!(st, "H2O@", wc * cmass / mtot * u"kg")
    V = volume(st)
    set_quantity!(st, "H+", 1e-7u"mol/L" * V.liquid)
    set_quantity!(st, "OH-", 1e-7u"mol/L" * V.liquid)
    return st
end

cement_mass(st) = sum(ustrip(us"kg", st.n[i] * cs.species[i][:M]) for i in ic)

function budget(α)
    fresh = paste(1.0)
    mc    = cement_mass(fresh)
    w_tot = ustrip(us"mol", fresh.n[iw]) * M_H2O
    eq, cert = equilibrate_certified(paste(α))

    n = collect(eq.n)                       # put the unreacted cement back
    for i in ic
        n[i] += (1 - α) * fresh.n[i]
    end
    ϕ = porosity(ChemicalState(cs, n), fresh)

    w_free = ustrip(us"mol", eq.n[iw]) * M_H2O
    V_ref  = ustrip(us"m^3", volume(fresh).total)
    return (; b = (w_tot - w_free) / mc / α,
              s = ϕ.void * V_ref / mc / α * 1e3,      # m³/kg → cm³/g
              w_free = w_free / mc, porosity = ϕ.total, certified = cert.optimal)
end

println(" α      b (g/g)   s (cm³/g)   free water (g/g)   porosity   certified")
for α in (0.55, 0.60, 0.65, 0.70, 0.80)
    r = budget(α)
    @printf("%5.2f  %9.4f  %10.4f  %16.4f  %9.4f   %s\n",
            α, r.b, r.s, r.w_free, r.porosity, r.certified)
end
```

Two things to read off that table.

**``b`` and ``s`` are constant in ``\alpha``** to four digits. That is not
assumed anywhere — it is what makes the linear budget of §1 legitimate for this
system, and it is checked by dividing each by ``\alpha`` and seeing the same
number five times.

**Every solve is certified.** For a convex problem the KKT conditions are
sufficient, so `certified = true` is a proof that the composition is the Gibbs
minimum and not the point an iteration stopped at — see
[Proving that an answer is the answer](@ref sec-theory-certificate).

```@example sd
ref = budget(0.65)
b_model, s_shrink = ref.b, ref.s
@printf("b = %.4f g/g (formula water)     s = %.4f cm³/g (chemical shrinkage)\n",
        b_model, s_shrink)
```

For scale: a chemical shrinkage of 0.064 cm³ per gram of cement is the value
cement chemistry reports for an ordinary Portland cement, and it is here a
**consequence** of the standard molar volumes in the database rather than an
input.

---

## 3. The measurement: ``S^\ast`` from a desorption isotherm

The retention curve is not chemistry. It says how tightly a particular material
holds the water still in it, it is measured, and it is the one thing on this page
that cannot come out of a thermodynamic database.

[BaroghelBouny1999](@cite) measured water-vapor desorption isotherms on two
pastes and two concretes and fitted each with

```math
p_c(S) = a\left(S^{-b} - 1\right)^{1 - 1/b},
```

which is the van Genuchten form written with ``b = 1/m`` — they say so where they
introduce Mualem's relative permeability. So a ``b`` from that literature enters
[`VanGenuchten`](@ref) as ``m = 1/b``, and **getting that inversion wrong is
silent**: it changes the exponent and nothing complains.

Their Table 5, with the mixes from their Table 1:

| mix | material | W/C | ``a`` (MPa) | ``b`` | ``m = 1/b`` |
|:--|:--|:--|:--|:--|:--|
| CO | cement paste | 0.34 | 37.5479 | 2.1684 | 0.46117 |
| CH | paste, 10 % silica fume | 0.19 | 96.2837 | 1.9540 | 0.51177 |
| BO | concrete | 0.48 | 18.6237 | 2.2748 | 0.43960 |
| BH | concrete, 10 % silica fume | 0.26 | 46.9364 | 2.0601 | 0.48541 |

The capillary pressure becomes a water activity through Kelvin,
``a_w = \exp(-p_c V_m / RT)``, which is what [`water_activity`](@ref) does:

```@example sd
co    = VanGenuchten(; a = 37.5479e6, m = 1 / 2.1684)   # their mix CO
V_m   = 1.807e-5      # m³/mol, liquid water at 25 °C
T_K   = 298.15
γ_w   = 0.0728        # N/m

a_w(S) = water_activity(co, S; V_m = V_m, T = T_K)

function saturation_at(rh)          # the form is not invertible in closed form
    lo, hi = 1e-4, 1 - 1e-12
    for _ in 1:60
        mid = (lo + hi) / 2
        a_w(mid) > rh ? (hi = mid) : (lo = mid)
    end
    return (lo + hi) / 2
end

println("   RH      S*     p_c (MPa)   Kelvin radius (nm)")
for rh in (0.75, 0.80, 0.85, 0.90, 0.95)
    S = saturation_at(rh)
    @printf("  %4.2f  %6.4f  %10.2f  %18.2f\n", rh, S,
            capillary_pressure(co, S) / 1e6,
            kelvin_radius(rh; γ = γ_w, V_m = V_m, T = T_K) * 1e9)
end
```

The last column sets the scale: 80 % relative humidity corresponds to a meniscus of
radius 4.8 nm. That is the **gel-pore scale**, which is why the water Powers
assigns to gel pores and the water a sealed paste cannot use are the same water.

```@example sd
using Plots

Ss = range(0.40, 0.999; length = 200)
p1 = plot(Ss, a_w.(Ss); xlabel = "degree of saturation S", ylabel = "water activity a_w",
    label = "mix CO, measured fit", linewidth = 2, color = :steelblue,
    title = "Desorption isotherm (Baroghel-Bouny et al. 1999)", legend = :bottomright)
hline!(p1, [0.80]; label = "RH = 0.80", linestyle = :dash, color = :firebrick)
vline!(p1, [saturation_at(0.80)]; label = "S* = $(round(saturation_at(0.80); digits = 4))",
    linestyle = :dash, color = :seagreen)
plot(p1; size = (700, 420), left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)
```

---

## 4. Closing the budget

Everything the boxed formula needs is now in hand, each from its own source.

```@example sd
powers_k(b, s, S★) = b + s * S★ / (1 - S★)

S80 = saturation_at(0.80)
@printf("S* = %.4f  so  S*/(1-S*) = %.4f\n\n", S80, S80 / (1 - S80))
@printf("with the model's formula water, b = %.4f:\n", b_model)
@printf("   k = %.4f      α_max(w/c=0.34) = %.4f     (Powers: 0.42 and %.4f)\n",
        powers_k(b_model, s_shrink, S80), wc / powers_k(b_model, s_shrink, S80),
        min(1.0, wc / 0.42))
```

**30 % above Powers.** Before calling that a disagreement, look at what ``b``
means on each side.

Powers' 0.42 is itself a sum: about 0.23 g/g of **non-evaporable** water plus
about 0.19 g/g of **gel** water. And his 0.23 is an *operational* quantity — the
water that survives D-drying, over ``\ce{Mg(ClO4)2*2H2O}`` at roughly 8 %
relative humidity. A thermodynamic database draws the line elsewhere: CEMDATA18
writes interlayer water **into the C-S-H formula**, where D-drying would have
removed it. The model's ``b`` and Powers' ``w_n`` are therefore not the same
quantity, and their difference is not an error:

```@example sd
b_powers = 0.23                     # Powers' non-evaporable water, his own number
@printf("model formula water   %.4f g/g\n", b_model)
@printf("Powers' w_n           %.4f g/g   (defined by D-drying)\n", b_powers)
@printf("difference            %.4f g/g   — interlayer water, counted differently\n",
        b_model - b_powers)
```

Substituting Powers' own non-evaporable water, and leaving the chemical shrinkage
and the isotherm exactly as measured:

```@example sd
@printf("with Powers' w_n = %.2f:\n", b_powers)
@printf("   k = %.4f      α_max(w/c=0.34) = %.4f     (Powers: 0.42 and %.4f)\n",
        powers_k(b_powers, s_shrink, S80), wc / powers_k(b_powers, s_shrink, S80),
        min(1.0, wc / 0.42))
```

11 % above. The remaining gap is a fifth of a point of saturation, and §6 shows
how little that is.

---

## 5. The result: inverting Powers

The comparison above runs one way — assume the arrest is at 80 % and see what
``k`` comes out. Run it the other way and the answer is a *prediction of a
measurable quantity*: **at what internal humidity does Powers' 0.42 place the
arrest?**

```@example sd
function invert_k(b, s, k_target)
    lo, hi = 1e-6, 1 - 1e-12
    for _ in 1:80
        mid = (lo + hi) / 2
        powers_k(b, s, mid) > k_target ? (hi = mid) : (lo = mid)
    end
    return (lo + hi) / 2
end

println("k = 0.42 implies:")
for (label, b) in (("the model's formula water", b_model), ("Powers' own w_n", b_powers))
    S = invert_k(b, s_shrink, 0.42)
    @printf("   b = %.4f (%-26s)   S* = %.4f   RH = %.4f\n", b, label, S, a_w(S))
end
```

With Powers' own non-evaporable water, **Powers' 0.42 corresponds to a hydration
arrest at 77.5 % relative humidity.**

That number was not put in anywhere. It came out of a chemical shrinkage computed
from a thermodynamic database, a desorption isotherm measured for a drying study,
and Powers' own water split — three sources, none of which knows what the others
are for. And 75–80 % is the window in which sealed pastes are independently
reported to stop hydrating; it is the same window `humidity_factor` implements as
its cut-off.

!!! tip "A corroboration from the same paper's own measurements"
    [BaroghelBouny1999](@cite) also measured the internal relative humidity of
    their sealed specimens at 28 days, in their Table 3: **97 %** for the paste
    CO (W/C 0.34) and **88.5 %** for the paste CH (W/C 0.19); **97 %** for the
    concrete BO (W/C 0.48) and **77.5 %** for the concrete BH (W/C 0.26).

    Monotone in W/C within each material class, as self-desiccation requires. And
    the mixes Powers says should have arrested by 28 days — the low-W/C ones —
    are the ones sitting in the 77.5–88.5 % band that the inversion above points
    at, while the two mixes with water to spare sit at 97 %. Those specimens are
    not arrested pastes at equilibrium, so this is corroboration of a range and
    not a fourth decimal; it is worth stating because the numbers come from the
    same table as the isotherm.

---

## 6. How much of this is the isotherm?

``k`` depends on ``S^\ast``, and ``S^\ast`` depends on the assumed arrest
humidity. The derivative is worth knowing before believing any single number:

```math
\frac{\mathrm{d}k}{\mathrm{d}S^\ast} = \frac{s}{(1 - S^\ast)^2}.
```

```@example sd
@printf("dk/dS* = %.3f at S* = %.4f\n\n", s_shrink / (1 - S80)^2, S80)
println("   RH      S*        k (model b)   k (Powers w_n)   α_max(0.34)")
for rh in (0.70, 0.75, 0.80, 0.85, 0.90)
    S = saturation_at(rh)
    @printf("  %4.2f  %6.4f  %13.4f  %16.4f  %12.4f\n", rh, S,
            powers_k(b_model, s_shrink, S), powers_k(b_powers, s_shrink, S),
            min(1.0, wc / powers_k(b_powers, s_shrink, S)))
end
```

```@example sd
rhs = range(0.65, 0.92; length = 120)
ks_m = [powers_k(b_model, s_shrink, saturation_at(r)) for r in rhs]
ks_p = [powers_k(b_powers, s_shrink, saturation_at(r)) for r in rhs]

p2 = plot(rhs, ks_m; xlabel = "assumed arrest humidity", ylabel = "k = (w/c) / α_max",
    label = "b = model formula water", linewidth = 2, color = :steelblue,
    title = "What the coefficient depends on", legend = :topleft, ylims = (0.25, 1.0))
plot!(p2, rhs, ks_p; label = "b = Powers' w_n = 0.23", linewidth = 2, color = :firebrick)
hline!(p2, [0.42]; label = "Powers 0.42", linestyle = :dash, color = :black)
hline!(p2, [0.36]; label = "Powers 0.36, saturated curing", linestyle = :dot, color = :gray)
plot(p2; size = (700, 420), left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)
```

A tenth of a point of saturation moves ``k`` by 0.14, so the arrest humidity has
to be known to about one point to pin ``k`` to 0.01. This page cannot separate
Powers' sealed 0.42 from his saturated-curing 0.36 on the strength of the
isotherm alone, and it does not try.

---

## 7. Two negative controls

Two claims carry the rest of the page: that the arrest criterion belongs in the
rate law and not in the Gibbs energy, and that Powers' proportional *form* is no
evidence for his coefficient. Both are testable on the material at hand, so both
are measured here rather than asserted.

### The thermodynamic route, measured

[`CapillaryWater`](@ref) imposes a water activity on the equilibrium itself. If
self-desiccation arrested hydration thermodynamically, imposing the arrest
humidity on a paste with all its cement available would leave some of that cement
unreacted. It does not:

```@example sd
fresh = paste(1.0)
ip = findfirst(s -> symbol(s) == "Portlandite", cs.species)
ij = findfirst(s -> symbol(s) == "Jennite", cs.species)

println("  a_w imposed   n(Portlandite)   n(Jennite)   free water (mol)   certified")
for aw in (1.00, 0.90, 0.80, 0.50)
    eq, cert = equilibrate_certified(
        paste(1.0); constraint = CapillaryWater(_ -> aw; reference = fresh)
    )
    @printf("  %11.2f   %14.6f   %10.6f   %16.5f   %s\n", aw,
            ustrip(us"mol", eq.n[ip]), ustrip(us"mol", eq.n[ij]),
            ustrip(us"mol", eq.n[iw]), cert.optimal)
end
```

From saturation down to the arrest humidity the assemblage is identical to six
digits — the imposed shift moves the answer by less than the printing precision —
and each of those answers is certified. The last row is there to show the other
half of that: further down, no route certifies, the warning above the table is
[`equilibrate_certified`](@ref) saying so, and the answer it returns is a KKT
candidate rather than a proof. The arithmetic in the opening admonition says why
nothing is expected to happen in that range anyway.

So a Gibbs minimization under [`CapillaryWater`](@ref) arrests where it ran out
of water stoichiometrically, not where Powers says. What the constraint is good
for is what §3–§5 use it for: making the water activity of a state mean the water
in the pores rather than a mole fraction, and handing that humidity to a rate law
through [`PoreHumidity`](@ref).

### The proportional form, measured

§1 argued that ``\alpha_{\max} \propto w/c`` follows from a fixed humidity
threshold composed with *any* monotone retention curve. The test is to run the
same construction on a second measured curve from the same table — mix BH, a
different material at a different W/C — and watch which number moves.

```@example sd
bh = VanGenuchten(; a = 46.9364e6, m = 1 / 2.0601)   # their mix BH

function saturation_of(law, rh)
    lo, hi = 1e-4, 1 - 1e-12
    for _ in 1:60
        mid = (lo + hi) / 2
        water_activity(law, mid; V_m = V_m, T = T_K) > rh ? (hi = mid) : (lo = mid)
    end
    return (lo + hi) / 2
end

for (name, law) in (("CO", co), ("BH", bh))
    S = saturation_of(law, 0.80)
    k = powers_k(b_model, s_shrink, S)
    @printf("\ncurve %s:  S* = %.4f   k = %.4f\n", name, S, k)
    for w in (0.25, 0.30, 0.35, 0.40)
        @printf("   w/c %.2f  →  α_max = %.6f   α_max/(w/c) = %.6f\n", w, w / k, 1 / k)
    end
end
```

The last column is the same number at every ``w/c``, in both blocks, to every
digit printed — and the two blocks disagree with each other by 18 %. That is the
whole point: the *constancy* of ``\alpha_{\max}/(w/c)`` is an identity of the
construction and carries no information, while its *value*, ``1/k``, carries all
of it. A page that quoted a five-digit agreement between that ratio and
``1/0.42 = 2.381`` would be quoting its own arithmetic, which is why §5 reports
``k`` and §6 reports ``\mathrm{d}k/\mathrm{d}S^\ast`` instead.

Neither curve gives 0.42 at RH 0.80 with this system's formula water. Recovering
Powers' coefficient takes his own water split as well, and §5 does that
explicitly and says so.

---

## 8. What went in, and what came out

| input | value | source | is that source about Powers? |
|:--|:--|:--|:--|
| clinker composition | Table 2 | [BaroghelBouny1999](@cite) | no |
| thermodynamic data | Cemdata18 | [Lothenbach2019](@cite) | no |
| ``b``, ``s`` | 0.3095 g/g, 0.0639 cm³/g | computed here, certified | no |
| retention curve | ``a`` = 37.5479 MPa, ``b`` = 2.1684 | [BaroghelBouny1999](@cite) Table 5 | no |
| ``\gamma``, ``V_m``, ``T`` | 0.0728 N/m, 1.807e-5 m³/mol, 298.15 K | water at 25 °C | no |
| ``w_n`` = 0.23 | Powers' own split of his 0.42 | [Powers1948](@cite) | **yes** |

**output**: the arrest humidity, 77.5 %.

Only the last row knows about Powers, and it contributes his *decomposition*, not
his coefficient. Nothing on this page was adjusted to improve the agreement, and
the one number that could have been — the retention curve — was fitted by its
authors to a drying experiment two decades before this calculation existed.

---

## 9. Assumptions, and where each one bites

  - **The degree of hydration is imposed, not predicted.** §2's construction
    reacts a fraction of the cement and leaves the rest inert. What this page
    computes is the ``\alpha`` at which the arrest criterion is met, not a
    trajectory in time. For that, hand [`PoreHumidity`](@ref) to
    [`parrot_killoh_avrami`](@ref) as its `humidity` and integrate.
  - **One retention curve for an evolving pore structure.** The measured isotherm
    belongs to a mature paste; the page applies it at every ``\alpha``. The pore
    structure of a young paste is coarser, so its true ``S^\ast`` is higher and
    the arrest earlier.
  - **The published fit excludes the humidity range this page works in.** The
    authors state that "the experimental data corresponding to the lowest
    capillary pressures (corresponding to the highest RH) have not been accounted
    [for] due to their weak reliability". Above about 90 % RH the curve used here
    is extrapolation into a region its authors declined to fit — which is exactly
    where an unarrested paste sits, and a reason the corroboration in §5 is
    stated as a range.
  - **The saturations are normalized differently.** The paper's ``S`` is relative
    to its measured total porosity, 30.3 % for mix CO; this model gives 27.4 % at
    ``\alpha = 0.65`` and reaches 30.3 % nearer ``\alpha \approx 0.57``. A 10 %
    difference in the denominator is a real bias on ``S^\ast``, and by §6 worth
    a few hundredths of ``k``.
  - **The C-S-H model matters as much as the isotherm.** ``\mathrm{d}k/\mathrm{d}b = 1``,
    so a change in the C-S-H water content moves ``k`` one for one. Measured on
    this system: `Jennite` gives ``b`` = 0.3095 g/g and the `CSHQ` solid solution
    gives **0.3684**, so it moves *away* from Powers, not toward him. Neither is
    wrong; they draw the formula-water line in different places.
  - **No alkalis.** This species list carries none, so the pore solution's
    osmotic depression of ``a_w`` is absent. At 0.1–0.5 mol/kg it is worth one to
    two points of relative humidity, in the same direction as the capillary term.
  - **Ideal molar volumes**, and a closed species list — the standing assumptions
    of every equilibrium page here.

## 10. Reproducing this

```bash
julia --project=scripts scripts/self_desiccation_powers.jl
```

Every block above is executed when this page is built, so the numbers printed are
whatever the code produced. The script computes the same quantities outside
Documenter and prints the intermediate ones this page summarizes.

See also: [`CapillaryWater`](@ref), [`PoreHumidity`](@ref),
[`WaterRetention`](@ref), [`powers_alpha_max`](@ref), and
[the w/c example](@ref sec-wc-ratio), which shows what a Gibbs minimization does
*without* an arrest criterion.
