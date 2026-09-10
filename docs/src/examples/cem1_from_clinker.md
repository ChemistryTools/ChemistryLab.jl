# [A CEM I from its clinker phases](@id sec-cem1-from-clinker)

This page starts where a cement data sheet stops: four anhydrous phases and a
sulfate carrier, a water/cement ratio, and nothing else. Everything after that —
how fast each phase reacts, which hydrates appear and in what amounts, how the
pore space closes, how the paste dries itself, how much heat it releases — is
computed.

It is meant to be read as the entry point to the cement material in this
documentation. [The w/c example](@ref sec-wc-ratio) asks what the *equilibrium*
assemblage is at a given degree of hydration;
[Self-desiccation](@ref sec-self-desiccation) asks where hydration stops. Here
the question is simpler and comes first: what does a paste do, over ninety days,
starting from its clinker.

## 1. Three clinkers, one controlled difference

The reference composition is the CEM I of [BaroghelBouny1999](@cite), Table 2.
The other two are **constructed** from it rather than measured: they trade alite
for belite at constant total silicate, which isolates the one effect this page
is about. They are not three commercial cements, and nothing here should be read
as a claim about any.

```@example cem1
using ChemistryLab
using OrdinaryDiffEq
using DynamicQuantities
using OrderedCollections
using Plots
using Printf

# mass fractions of the anhydrous cement. C3S + C2S is held at 0.8126.
CLINKERS = OrderedDict(
    "measured (Baroghel-Bouny CO)" =>
        (C3S = 0.5728, C2S = 0.2398, C3A = 0.0303, C4AF = 0.0759, Gp = 0.0439),
    "alite-rich (constructed)" =>
        (C3S = 0.7126, C2S = 0.1000, C3A = 0.0303, C4AF = 0.0759, Gp = 0.0439),
    "belite-rich (constructed)" =>
        (C3S = 0.4128, C2S = 0.3998, C3A = 0.0303, C4AF = 0.0759, Gp = 0.0439),
)
WC = 0.45          # enough water that the arrest is not the subject here

for (name, c) in CLINKERS
    @printf("  %-30s C3S %.4f  C2S %.4f  C3A %.4f  C4AF %.4f  Gp %.4f\n",
            name, c.C3S, c.C2S, c.C3A, c.C4AF, c.Gp)
end
```

## 2. The chemical system

The species list is the clinker, the sulfate carrier and the hydrates they can
form. Nothing else is admitted, and that closed list is an assumption worth
stating: a phase absent from it cannot appear however stable it would be.

```@example cem1
substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
selection = split(
    "C3S C2S C3A C4AF Gp Portlandite Jennite ettringite monosulphate12 " *
        "C3AH6 C3FH6 H2O@"
)
cs = ChemicalSystem(
    speciation(substances, selection; aggregate_state = [AS_AQUEOUS]), CEMDATA_PRIMARIES
)
sp(n) = cs[n]
println("$(length(cs.species)) species in the system")
```

## 3. The reactions, and the rate law each one carries

Four hydration reactions, each with a Parrot-Killoh rate in its Avrami form.
`α_max = 1.0` and the arrest, if any, is left to the humidity coupling of §7.

```@example cem1
law = VanGenuchten(; a = 37.5479e6, m = 1 / 2.1684)   # Baroghel-Bouny, mix CO

function build(compo; humidity = true, tend = 90 * 86400.0)
    st = ChemicalState(cs)
    for (name, frac) in pairs(compo)
        set_quantity!(st, string(name), frac * u"kg")
    end
    set_quantity!(st, "H2O@", WC * sum(values(compo)) * u"kg")
    h = humidity ? PoreHumidity(law, cs; reference = st) : nothing

    specs = (
        ("C3S", PK84_PARAMS_C3S,
         OrderedDict(sp("C3S") => 1.0, sp("H2O@") => 103 / 30),
         OrderedDict(sp("Jennite") => 1.0, sp("Portlandite") => 4 / 3)),
        ("C2S", PK84_PARAMS_C2S,
         OrderedDict(sp("C2S") => 1.0, sp("H2O@") => 73 / 30),
         OrderedDict(sp("Jennite") => 1.0, sp("Portlandite") => 1 / 3)),
        # No gypsum in this one, deliberately — see §9.
        ("C3A", PK84_PARAMS_C3A,
         OrderedDict(sp("C3A") => 1.0, sp("H2O@") => 6.0),
         OrderedDict(sp("C3AH6") => 1.0)),
        ("C4AF", PK84_PARAMS_C4AF,
         OrderedDict(sp("C4AF") => 1.0, sp("Portlandite") => 2.0, sp("H2O@") => 10.0),
         OrderedDict(sp("C3AH6") => 1.0, sp("C3FH6") => 1.0)),
    )
    rxns = AbstractReaction[]
    for (nm, pk, reac, prod) in specs
        rx = Reaction(reac, prod; symbol = nm)
        rx[:rate] = parrot_killoh_avrami(
            pk, nm; α_max = 1.0, blaine = 380.0u"m^2/kg", humidity = h
        )
        push!(rxns, rx)
    end
    kp = KineticsProblem(cs, rxns, st, (0.0, tend))
    ks = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-10)
    return kp, integrate(kp, ks), st, h
end

runs = OrderedDict(name => build(c) for (name, c) in CLINKERS)
println("integrated: ", join(keys(runs), ", "))
```

## 4. How fast each phase reacts

The four clinker phases have very different kinetics, and the shape of the early
heat curve of a cement is mostly the aluminate's doing while its strength is
mostly the alite's.

```@example cem1
kp, sol, st0, h = runs["measured (Baroghel-Bouny CO)"]
ts = exp10.(range(log10(3600.0), log10(90 * 86400.0); length = 120))
dh = degrees_of_hydration(sol, kp; times = ts)
days = ts ./ 86400

p1 = plot(; xscale = :log10, xlabel = "time (days)", ylabel = "degree of hydration α",
    title = "Clinker phases, measured composition", legend = :topleft, ylims = (0, 1))
for (nm, col) in (("C3S", :steelblue), ("C2S", :seagreen),
                  ("C3A", :firebrick), ("C4AF", :darkorange))
    plot!(p1, days, dh[nm]; label = nm, linewidth = 2, color = col)
end
plot(p1; size = (720, 430), left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)
```

Alite and the aluminate move together early — both near 0.36 at one day — and
the ferrite and belite lag well behind, at 0.14 and 0.10. By twenty-eight days
alite is at 0.75 and belite has reached only 0.39; at ninety days they are 0.83
and 0.50. Belite is not so much a slow phase as a **late** one, and that
separation is why a cement keeps gaining strength for months rather than days.

## 5. What changing the clinker does

```@example cem1
p2 = plot(; xscale = :log10, xlabel = "time (days)",
    ylabel = "silicate degree of hydration",
    title = "Alite/belite balance, at constant total silicate",
    legend = :topleft, ylims = (0, 1))
for ((name, (kpi, soli, _, _)), col) in zip(runs, (:black, :firebrick, :seagreen))
    d = degrees_of_hydration(soli, kpi; times = ts)
    c = CLINKERS[name]
    blend = (c.C3S .* d["C3S"] .+ c.C2S .* d["C2S"]) ./ (c.C3S + c.C2S)
    plot!(p2, days, blend; label = name, linewidth = 2, color = col)
end
plot(p2; size = (760, 430), left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)
```

Trading alite for belite at constant silicate content separates the three curves
by about 0.08 in degree of hydration at one day (0.244 to 0.328) — and that
absolute spread does **not** close: it is still 0.096 at ninety days (0.694 to
0.790). What closes is the spread *relative* to how far the reaction has gone,
from 29 % of the mean at one day to 13 % at ninety.

So the alite-rich cement is ahead at every age and is never overtaken. The
familiar statement that a belite cement "catches up" is a statement about the
relative gap, not the absolute one, and the two are worth keeping apart.

## 6. The hydrates that appear

```@example cem1
states = [state_at(sol, kp, t) for t in ts]
amount(name) = [ustrip(us"mol", s.n[findfirst(x -> symbol(x) == name, cs.species)])
                for s in states]

p3 = plot(; xscale = :log10, xlabel = "time (days)", ylabel = "amount (mol)",
    title = "Hydrate assemblage, measured composition", legend = :topleft)
for (nm, lbl, col) in (("Portlandite", "portlandite CH", :steelblue),
                       ("Jennite", "C-S-H (Jennite)", :seagreen),
                       ("C3AH6", "hydrogarnet C3AH6", :firebrick),
                       ("C3FH6", "ferrite hydrogarnet", :darkorange))
    plot!(p3, days, amount(nm); label = lbl, linewidth = 2, color = col)
end
plot(p3; size = (720, 430), left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)
```

The two silicate hydrates — portlandite and C-S-H — keep growing for the whole
ninety days, which is §4's alite and belite curves seen from the product side.
The aluminate hydrogarnet is nearly complete by a few days, following its own
clinker phase. **No ettringite appears at all**, and the gypsum in the mix is
never touched: §9 is about why.

## 7. The pore space closes, and the paste dries itself

Hydrates occupy less than the water and clinker they were made from, so voids
open even as the solids grow: that is chemical shrinkage, and it is what empties
the pores of a sealed paste.

```@example cem1
V0 = volume(st0)
poro = [porosity(s, st0).total for s in states]
# m³ → cm³ by hand: `us"cm^3"` is a symbolic unit and does not convert a
# quantity carried in SI dimensions.
shrink = [1.0e6 * ustrip(us"m^3", V0.total - volume(s).total) for s in states]
rh = [h(ustrip.(us"mol", s.n)) for s in states]

p4 = plot(days, 100 .* poro; xscale = :log10, xlabel = "time (days)",
    ylabel = "total porosity (%)", label = "porosity", linewidth = 2,
    color = :steelblue, title = "Pore space and self-desiccation", legend = :left)
p4b = twinx(p4)
plot!(p4b, days, rh; xscale = :log10, ylabel = "internal relative humidity",
    label = "internal RH", linewidth = 2, color = :firebrick, legend = :right,
    ylims = (0.75, 1.01))
hline!(p4b, [0.80]; label = "humidity_factor cut", linestyle = :dash, color = :gray)
plot(p4; size = (760, 430), left_margin = 10Plots.mm, right_margin = 12Plots.mm,
    bottom_margin = 8Plots.mm)
```

```@example cem1
p5 = plot(days, shrink; xscale = :log10, xlabel = "time (days)",
    ylabel = "volume lost (cm³)", label = "chemical shrinkage", linewidth = 2,
    color = :seagreen, title = "Chemical shrinkage", legend = :topleft)
plot(p5; size = (720, 400), left_margin = 10Plots.mm, bottom_margin = 8Plots.mm)
```

At `w/c = 0.45` the internal humidity falls but does not reach the cut, so this
paste is not arrested by self-desiccation within ninety days.
[Self-desiccation in time](@ref sec-self-desiccation-kinetics) runs the drier
mixes where it is.

## 8. Why the aluminate reaction here has no gypsum in it

A real cement's aluminate does not go to hydrogarnet: while sulfate lasts it
makes **ettringite**, and once the sulfate is exhausted the ettringite converts
to monosulfate. Writing that first stage as a fixed-stoichiometry kinetic
reaction is tempting,

```julia
C3A + 3 Gp + 26 H2O → ettringite          # do not do this with a PK rate
```

and it is wrong for a reason worth showing rather than asserting. A
Parrot-Killoh rate depends on the degree of reaction of **its own clinker
phase** and on nothing else. It does not watch the gypsum. Running exactly that
reaction on this cement — as a variant, not on this page, since it is the thing
to avoid — gives:

| quantity | value |
|:--|:--|
| gypsum present at ``t = 0`` | 0.25499 mol |
| C₃A consumed by 90 days | 0.09388 mol |
| gypsum the extent demands (``3\times``) | 0.28164 mol |
| gypsum actually consumed | 0.25499 mol, and it stops at exactly zero |

The extent kept advancing after the gypsum ran out. The amount of gypsum is
floored at zero rather than going negative, so the shortfall — 0.027 mol of
sulfate — is simply **created**. Nothing in the package objects:
[`extent_residual`](@ref) measures the drift between the integrated species and
the integrated extents, which is unaffected, and the feasibility machinery
guards the element balance of an *equilibrium* sub-solve, which is not running
here.

So a fixed-stoichiometry kinetic reaction is only safe when its co-reactants
cannot run out. Two ways round it, and this page takes the first:

  - **write the reaction without the limiting co-reactant**, as here and as
    [Cement clinker hydration kinetics](@ref) does. The aluminate then goes to
    hydrogarnet, the sulfate stays untouched, and the page does not pretend to
    model the AFt/AFm sequence;
  - **let the assemblage be a result**: attach an equilibrium solver, so the
    kinetics supplies element budgets and the products are whatever minimizes
    the Gibbs energy at each step, ettringite while sulfate lasts and monosulfate
    afterwards. That is [Coupled hydration](@ref), and it is the honest route to
    a sulfate-bearing aluminate.

## 9. What this page assumed

  - **a closed species list** — the hydrates admitted in §2 and no others;
  - **stoichiometric reactions written by hand**, one per clinker phase, rather
    than an equilibrium assemblage. [Coupled hydration](@ref) attaches an
    equilibrium solver instead, and then which hydrates appear becomes a result;
  - **Parrot-Killoh parameters** as published, with a Blaine correction and no
    fitting to any measurement on this page;
  - **two constructed clinkers**, which exist to isolate one variable and are
    not cements anyone has made.

See also: [Cement clinker hydration kinetics](@ref) for the same machinery on
one composition with calorimetry, [the w/c example](@ref sec-wc-ratio) for the
equilibrium view, and [the water budget](@ref sec-theory-water-budget) for what
the arrest is and is not.
