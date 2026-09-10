# [The water budget of a hydrating paste](@id sec-theory-water-budget)

Mix one gram of cement with ``w/c`` grams of water, seal it, and wait. Some of
the water ends up inside hydrate formulas, some stays in pores the hydrates
leave behind, and the reaction stops well before the cement is gone. Powers'
rule of thumb summarizes the arrest as

```math
\alpha_{\max} = \frac{w/c}{0.42} ,
```

and [`powers_alpha_max`](@ref) supplies it to the rate laws. This page is about
what that 0.42 is made of, which part of it a Gibbs energy minimization can
predict, and — the question worth being precise about — **why the other part is
not a thermodynamic quantity at all**.

## 1. Three destinations for the mixing water

At a degree of hydration ``\alpha``, the mixing water has gone to three places:

| destination | symbol | what it is |
|:--|:--|:--|
| **chemically bound** | ``b\,\alpha`` | water written into hydrate formulas — C-S-H, portlandite, AFm, AFt. It is no longer water |
| **the volume the reaction lost** | ``s\,\alpha`` | hydrates occupy less than the water and clinker they came from, so voids open: **chemical shrinkage** |
| **still liquid** | the remainder | pore solution, held more or less tightly depending on the pore it sits in |

Both coefficients are **computed**, not assumed. From certified equilibria on a
CEM I at imposed ``\alpha`` — the calculation is executed on the
[self-desiccation page](@ref sec-self-desiccation) — this species list gives

```math
b = 0.3095 \ \mathrm{g/g} , \qquad s = 0.0639 \ \mathrm{cm^3/g} ,
```

each **constant in ``\alpha`` to four digits**, which is what makes a budget
linear in ``\alpha`` legitimate for this system rather than merely convenient.
``s \approx 0.064`` cm³ per gram of cement is the chemical shrinkage that cement
chemistry reports for an ordinary Portland cement, and here it is a consequence
of the standard molar volumes in the database.

## 2. What the minimization predicts: the stoichiometric threshold

A Gibbs minimization under a conservation budget answers one question: what is
the lowest-energy composition reachable *with the atoms available*. Water
appears in it as hydrogen and oxygen, and nothing more. So the threshold it
predicts is the point at which the hydrogen runs out — that is, ``w/c \approx b``.

Measured on the [w/c example](@ref sec-wc-ratio)'s species list, at equilibrium:

| ``w/c`` | clinker left at equilibrium |
|:--|:--|
| 0.15 | 54.6 % |
| 0.28 | 2.0 % |
| 0.30 and above | 0 % |

That is a real prediction and it is a **stoichiometric** one. Below ``w/c \approx
0.30`` clinker survives because there is not enough water to convert it, and the
number 0.30 is nothing but ``b`` plus the little the solution retains.

!!! warning "Below that threshold the calculation also leaves its own domain"
    At ``w/c = 0.28`` the free water goes to the solver's floor, the solvent mole
    fraction falls to 0.21 and the ionic strength reaches hundreds of mol/kg. The
    composition is arithmetically a solution and physically not one — every
    molality is per kilogram of a solvent that is no longer there. This is what
    [`solvent_fraction`](@ref) and [`SOLVENT_FRACTION_FLOOR`](@ref) exist to
    name, and no activity model repairs it.

## 3. What Powers' 0.42 actually contains

Powers splits his coefficient himself: about **0.23 g/g** of non-evaporable water
— what survives D-drying, his operational definition of chemically bound — plus
about **0.19 g/g** of **gel water**, which is water *present in the paste and
unusable by the reaction*.

So the gap between the two thresholds,

```math
0.30 \;\;(\text{stoichiometric, predicted})
\qquad\longleftrightarrow\qquad
0.42 \;\;(\text{Powers, measured}) ,
```

**is exactly the gel water**, and it is the whole of the question. A paste at
``w/c = 0.35`` has enough hydrogen to consume all its clinker — the minimization
says so, correctly — and does not do it.

Written as a budget, the arrest condition is

```math
\frac{w}{c} = \underbrace{b\,\alpha_{\max}}_{\text{bound}}
            + \underbrace{s\,\alpha_{\max}\frac{S^\ast}{1-S^\ast}}_{\text{the water left at arrest}} ,
\qquad\text{so}\qquad
\alpha_{\max} = \frac{w/c}{k} ,
\quad k = b + s\,\frac{S^\ast}{1-S^\ast} ,
```

where ``S^\ast`` is the degree of saturation at which hydration stops. Two of the
three ingredients are thermodynamic and computed; ``S^\ast`` is neither.

## 4. Why the gel water is not a thermodynamic effect

The natural hypothesis is that confined water is thermodynamically unavailable:
water held in a fine pore has a reduced activity, so hydrates that consume it
become less stable, and the arrest should fall out of the minimization.
**The magnitude is wrong by two orders of magnitude**, and the arithmetic is
short enough to do here:

```@example wb
using Printf
R, T = 8.31446261815324, 298.15
a_w = 0.80                     # the internal humidity a sealed paste arrests at
ΔG_water = R * T * log(a_w)                    # J per mole of water
n_water_per_alite = 3.3                        # mol H₂O per mol C3S → C-S-H + CH
@printf("RT ln a_w                        = %8.1f J/mol of water\n", ΔG_water)
@printf("× %.1f mol water per mol alite    = %8.2f kJ per mol of alite\n",
        n_water_per_alite, n_water_per_alite * ΔG_water / 1000)
@printf("hydration ΔG, order of magnitude = %8.0f kJ per mol of alite\n", -100.0)
@printf("\nactivity that would null it      = %.2e\n", exp(-100_000 / (n_water_per_alite * R * T)))
```

An activity of that size corresponds, through Kelvin, to a meniscus radius
**smaller than a water molecule**. And the measurement agrees with the
arithmetic: imposing a water activity from saturation down to 0.80 through
[`CapillaryWater`](@ref) leaves the equilibrium assemblage of a CEM I paste
unchanged to six digits, with a certificate on every answer — that is the first
negative control on the [self-desiccation page](@ref sec-self-desiccation).

So the capillary lowering of the water activity is real, representable, and
energetically negligible. What stops a real paste is elsewhere.

### It is a question of scale, which a 0D calculation does not have

The reaction needs water *at the surface of a clinker grain*, and dissolved ions
carried *away* from it. As hydration proceeds:

  - the coarse capillary pores empty first, because they hold water most loosely;
  - what remains sits in gel pores — nanometers wide, inside the C-S-H — as films
    one or two molecules thick;
  - the permeability of the paste falls by orders of magnitude, and the liquid
    path between the remaining water and an unreacted grain, tens of microns
    away, is broken.

The reactants stop meeting. That is a statement about **geometry and transport
across four orders of magnitude of length**, and a 0D model has no length at
all: it assumes one well-mixed solution in contact with every solid, so a water
molecule in a 1 nm gel pore is the same object as one in a 10 μm capillary pore.
No activity model, however good, can distinguish them — the distinction is not
thermodynamic.

There is a numerical hint of that boundary inside the theory itself. The
[Debye screening length](@ref sec-theory-activity) of a cement pore solution at
``I \approx 0.3`` mol/kg is **0.55 nm** — the thickness of two water molecules,
and the width of the pores in question. Where the screening length and the pore
are the same size, the continuum-dielectric and mean-field assumptions that
produce every activity coefficient on that page are both strained. The formulas
keep returning numbers; they have simply left the picture they were derived in.

## 5. So what *can* be predicted, and how

Three routes, and they differ in what has to be assumed.

**Impose ``\alpha``.** Give the minimization a degree of hydration and let it
compute the assemblage — the construction of [LothenbachWinnefeld2006](@cite),
used by the [w/c example](@ref sec-wc-ratio) below its stoichiometric bound and
by the self-desiccation page throughout. Honest, and ``\alpha`` is an input:
either measured, or taken from `powers_alpha_max`. This is also what GEM-Selektor
and Reaktoro offer, since they are 0D equilibrium codes too.

**Close the budget with a measured isotherm.** Take ``b`` and ``s`` from the
minimization and ``S^\ast`` from a published desorption isotherm at an assumed
arrest humidity. Then ``k`` is predicted up to that one empirical threshold — and
inverting it is more informative than quoting it: with Powers' own
``w_n = 0.23``, ``k = 0.42`` corresponds to an arrest at **77.5 % relative
humidity**, the window sealed pastes are independently reported to stop in. The
[self-desiccation page](@ref sec-self-desiccation) does this, and shows why the
*proportionality* ``\alpha_{\max}\propto w/c`` is structural and therefore no
evidence at all.

**Integrate a rate law that reads the humidity.** [`humidity_factor`](@ref)
implements the empirical cut that stops hydration below about 80 % RH, and
[`PoreHumidity`](@ref) computes that humidity from the current saturation, so a
kinetic run arrests on its own. The arrest becomes a result of the integration
rather than a criterion applied afterwards — while the *criterion* stays
empirical. That is the honest ceiling of a 0D framework.

**Beyond 0D.** Making the arrest a genuine prediction needs a model with a length
in it: a microstructural one, where the arrest emerges because the liquid path
percolates no longer, or a transport calculation on the pore network. That is a
different package.

## See also

  - [Self-desiccation](@ref sec-self-desiccation) — the budget closed, with the
    numbers computed and the negative controls executed
  - [The w/c ratio](@ref sec-wc-ratio) — the stoichiometric threshold, scanned
  - [Activity models](@ref sec-theory-activity) — why the screening length
    matters here
  - [`powers_alpha_max`](@ref), [`PoreHumidity`](@ref), [`CapillaryWater`](@ref)
