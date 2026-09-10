# [Activity models](@id sec-theory-activity)

Every equilibrium this package computes rests on

```math
\mu_i = \mu_i^\circ + RT \ln a_i ,
```

and on nothing else about the solution. The standard potential ``\mu_i^\circ``
comes from the thermodynamic database; the activity ``a_i`` comes from the
**activity model**, which is therefore the single place where a real solution
stops being ideal. This page is about what those models are, where they come
from, what they need, and where they stop being true.

A model has to answer two questions, not one:

1. the **solute** activity ``a_i = \gamma_i m_i``, which sets every solubility
   and every saturation index;
2. the **solvent** activity ``a_w``, which is what a hydrate reaction consumes,
   and which in a cement paste short of mixing water is the quantity that
   decides how far the reaction can go.

The three built-in models — [`DiluteSolutionModel`](@ref),
[`DaviesActivityModel`](@ref), [`HKFActivityModel`](@ref) — differ in both, and
§5 measures by how much.

## 1. Where the ``\sqrt{I}`` comes from

An ion in a solution of ions is not in the same state as one alone at the same
molality, and the reason is electrostatic. Around a cation, anions are slightly
more likely to be found than cations — not because of any structure, only
because their energy there is lower — so the ion sits at the center of a diffuse
**ionic cloud** of opposite net charge. That cloud screens the ion's own field
and lowers its energy, which is why activity coefficients of ions are *below* 1.

The classical treatment takes three steps, and the assumptions it makes on the
way are exactly the assumptions that later fail:

1. **Poisson's equation** relates the mean potential to the mean charge density,
   ``\nabla^2\psi = -\rho/(\varepsilon_0\varepsilon_r)``, with the solvent
   entering only as a **continuum of bulk permittivity** ``\varepsilon_r``.
2. **Boltzmann statistics** give the local ion densities,
   ``n_j(r) = n_j^0\exp(-z_j e\psi/k_BT)``, treating the ions as *independent* in
   the mean field of all the others.
3. **Linearization**, ``z_j e \psi \ll k_B T``, turns the pair into
   ``\nabla^2\psi = \kappa^2\psi``, whose solution is a screened Coulomb
   potential with a single length scale

```math
\kappa^{-1} = \left(\frac{\varepsilon_0\varepsilon_r k_B T}
                        {2 e^2 N_A \rho_w I}\right)^{1/2}
```

the **Debye length**. Charging the central ion reversibly inside its own cloud
gives the excess chemical potential, and hence

```math
\log_{10}\gamma_i = -A z_i^2 \sqrt{I} \qquad\text{(the limiting law)} ,
```

exact as ``I \to 0``. The ``\sqrt{I}`` is not fitted: it is ``\kappa``, and
``\kappa \propto \sqrt{I}``.

Two ions cannot approach closer than the sum of their radii, so the cloud is
excluded from a shell of radius ``\mathring{a}`` around the central ion. Carrying
that through gives the **extended** law, which is what
[`HKFActivityModel`](@ref) implements:

```math
\log_{10}\gamma_i = -\frac{A z_i^2\sqrt{I}}{1 + B\mathring{a}_i\sqrt{I}} ,
\qquad
I = \tfrac{1}{2}\sum_j m_j z_j^2 .
```

### ``A`` and ``B`` are properties of water, not fitting constants

Collecting the constants of step 3 gives, with ``\rho_w`` in g/cm³ and
``\varepsilon`` the dielectric constant of water,

```math
A = 1.824829238\times10^{6}\,\frac{\sqrt{\rho_w}}{(\varepsilon T)^{3/2}} ,
\qquad
B = 50.29158649\,\frac{\sqrt{\rho_w}}{\sqrt{\varepsilon T}} ,
```

which is what [`hkf_debye_huckel_params`](@ref) evaluates from this package's own
equation of state for water. So the ``A = 0.5114`` and ``B = 0.3288`` that the
models carry as defaults are **derived**, not adopted, and they agree with
[Helgeson1981](@cite) Table 1 — evaluated at three temperatures in
[What the choice of activity model costs](@ref sec-app-activity-models).

Both rise with temperature, because water's dielectric constant falls faster than
``T`` rises: hot water screens worse, so the same ionic strength costs more.

### The screening length is the size of a gel pore

The Debye length is ``\kappa^{-1} = 1/(B\sqrt{I})`` in ångström when ``B`` is in
Å⁻¹(kg/mol)^½. At 25 °C it runs from
9.6 nm at ``I = 10^{-3}`` mol/kg to 0.30 nm at 1 mol/kg, passing **0.55 nm at
``I = 0.3`` mol/kg** — the values are tabulated in
[What the choice of activity model costs](@ref sec-app-activity-models).

A cement pore solution sits around ``I \approx 0.1``–``0.5 mol/kg``, so its
screening length is a **few ångström** — the thickness of two or three water
molecules. That is the same scale as the water films in the gel pores of C-S-H,
and both assumptions of step 1 and step 2 above are strained there: a continuum
of bulk permittivity, and ions independent in a mean field. Nothing in the
formulas announces it. It is the physical reason to treat an extended
Debye-Hückel model as a correlation valid in bulk solution rather than as a
theory of confined water.

## 2. The ``\dot{B} I`` term is a deviation function, not a physical term

Measured activity coefficients turn back *upwards* at high ionic strength, which
no screening argument produces. The B-dot model adds a linear term for it,

```math
\log_{10}\gamma_i = -\frac{A z_i^2\sqrt{I}}{1 + B\mathring{a}_i\sqrt{I}}
                    + \dot{B} I .
```

Its status is that of an empirical correlation. [AndersonCrerar1993](@cite) (§17.7.1,
pp. 445–446) record that Helgeson defined ``\dot{B}`` as a **deviation
function**: the difference between the *observed* activity coefficient of an
electrolyte — NaCl — and what the extended Debye-Hückel expression predicts for
it. So it carries short-range ion-solvent and ion-ion interaction *and* whatever
the first two terms failed to capture, together, in one number fitted to one
salt. They add that [Helgeson1981](@cite) later split it into a hydration term
from the Born equation and a residual short-range term.

That is why this model has a *ceiling* rather than an asymptote, and why the
ceiling is quoted vaguely as "about a molal": the term is not wrong so much as
it is standing in for physics it does not contain.

Neutral species get the **Setschenow** form, ``\log_{10}\gamma_i = K_n I``:
water engaged in the solvation shells of ions is water unavailable to solvate a
neutral molecule, so its activity rises with ionic strength and its solubility
falls. This is salting out, and `CO₂(aq)` is the case that matters for
carbonation.

## 3. The water activity, and why it cannot be assumed separately

The solvent is not a solute and its activity is not obtained by the same
formula. Two routes exist in this package.

**Raoult** — ``a_w = x_w``, the mole fraction — is what
[`DiluteSolutionModel`](@ref) and [`DaviesActivityModel`](@ref) use. It counts
molecules and knows nothing about what they are.

**The osmotic coefficient** ``\varphi`` is what [`HKFActivityModel`](@ref) uses:

```math
\ln a_w = -M_w \varphi \sum_j m_j ,
```

with ``\varphi`` obtained by integrating the Gibbs-Duhem relation over the same
``A``, ``B`` and ``\dot{B}`` that produced the ``\gamma_i``. That is the whole
point of it. At constant ``T`` and ``P``,

```math
\sum_i n_i \,\mathrm{d}\mu_i = 0 ,
```

which is not an optional refinement: it is the statement that the solvent and
the solutes are parts of one thermodynamic system. A model that corrects its
solutes and leaves its solvent ideal violates it by construction, and §5
measures by how much.

The one approximation in the B-dot route is that ``\varphi`` uses a single
charge-weighted mean radius,
``\mathring{a}_{\text{eff}} = \sum_i m_i z_i^2\mathring{a}_i / \sum_i m_i z_i^2``,
where the ``\gamma_i`` use per-ion radii. §5 measures that too.

## 4. Inputs and outputs, model by model

What each model **needs** and what it **returns** — the full parameter tables,
with the provenance of every default, are in the docstrings
([`HKFActivityModel`](@ref), [`DaviesActivityModel`](@ref)); this is the summary
that lets you choose.

| | [`DiluteSolutionModel`](@ref) | [`DaviesActivityModel`](@ref) | [`HKFActivityModel`](@ref) | [`PitzerActivityModel`](@ref) |
|:--|:--|:--|:--|:--|
| solute scale | molarity | molality | molality | molality |
| ``\gamma_i`` | ``\equiv 1`` | Davies | extended D-H + ``\dot{B} I`` | virial expansion |
| ``a_w`` | Raoult | Raoult | osmotic coefficient | osmotic coefficient |
| per-species data | none | none | ion radii ``\mathring{a}_i`` (tabulated, overridable) | **a parameter per ion pair and per triplet — caller input** |
| scalar inputs | none | ``A``, ``b``, ``b_n`` | ``A``, ``B``, ``\dot{B}``, ``K_n``, ``\mathring{a}_{\text{default}}`` | the shape constants ``\alpha_1``, ``\alpha_2``, ``b`` |
| ``T``, ``P`` dependence | none | ``A(T,P)`` on request | ``A(T,P)``, ``B(T,P)`` on request | ``A_\varphi(T,P)`` only; the ``\beta`` set is fitted at one temperature |
| returns | ``\ln a_i`` for every species | same | same | same |
| ``\gamma`` useful to | ``I \lesssim 0.01`` | ``I \lesssim 0.5`` | ``I \lesssim 1`` | the range its set was fitted over, a few mol/kg |
| Gibbs-Duhem consistent | approximately | **no** (§5) | to ``10^{-5}`` (§5) | **exactly, by construction** (§6) |

All three return the same object — a vector of ``\ln a_i`` indexed like
`cs.species`, covering solutes, solvent, pure crystals (``0``), gases and
solid-solution end-members — so they are interchangeable at every call site, and
`concentration_scale` tells the accessors which convention was used.

## 5. What the difference is worth, measured

The models are compared, on an imposed NaCl composition and without solving
anything, in
[What the choice of activity model costs](@ref sec-app-activity-models). Three
results from that page belong here, because they are about the theory rather
than about the numbers:

**The ideal model is already several percent off at a millimolal.** Ideality is
not a safe default that degrades gracefully; it is exact only in a limit.

**Davies and the B-dot model part company around a tenth molal**, and by
3 mol/kg Davies returns ``\gamma > 1`` while the B-dot model is still below 1.
That is the ``bI`` term of §2 taking over from the screening term — the ceiling
of a deviation function arriving, visible in the numbers.

**The values of ``a_w`` barely separate at all** — Raoult and the osmotic route
differ by a few parts in a thousand even at 3 mol/kg — while their
**derivatives** differ by four orders of magnitude. Measured as the Gibbs-Duhem
residual ``\lvert\sum_i n_i\,\mathrm{d}\mu_i\rvert`` along a dissolution at
1 mol/kg: ``1.9\times10^{-1}`` for Davies against ``6.7\times10^{-5}`` for the
B-dot model, with the ideal model at ``2.9\times10^{-2}`` in between.

The third result follows from §3: **Davies is less thermodynamically consistent
than assuming ideality.** Correcting the
solutes while leaving the solvent at ``a_w = x_w`` sets the two halves of one
model against each other, and a model can be *more* wrong for being *partly*
corrected. Since equilibrium is set by derivatives and not by values, a
disagreement invisible in ``a_w`` is decisive in ``\mu_w``.

The B-dot model's own defect — the single charge-weighted mean radius in its
osmotic coefficient — appears only along a composition change that holds ``I``
and ``\sum m`` fixed, where it measures a few parts in a thousand. That is the
direction `test/activities.jl` uses, which is why its tolerance is `5e-3` rather
than solver tolerance.

## 6. The ion-interaction model: a different kind of object

Everything above is a **corrected Debye-Hückel law**: one screening term derived
from electrostatics, one size correction, and one empirical term standing in for
the rest. [`PitzerActivityModel`](@ref) is not that. It starts from the excess
Gibbs energy and expands it as a **virial series in the molalities**, exactly as
the pressure of a non-ideal gas is expanded in its density —
[AndersonCrerar1993](@cite) (§17.8) derive it as a cluster expansion with
osmotic pressure in place of pressure, which is where the analogy is exact:

> the first term in a virial equation always represents ideal behavior; in the
> second term ``B_2`` represents the non-ideal contribution from pairwise
> interactions of molecules; ``B_3`` gives the interactions of triples.

Hence the shape of the parameter set: one coefficient per **ion pair**, one per
**triplet**, and no per-species radius at all.

```math
\frac{G^{\text{ex}}}{RT n_w} =
  f(I)
  + \sum_{c}\sum_{a} m_c m_a\, B_{ca}(I)
  + \sum_{c<c'} m_c m_{c'}\,\theta_{cc'}
  + \sum_{a<a'} m_a m_{a'}\,\theta_{aa'}
  + \text{triplets } (\psi)
  + \text{neutrals } (\lambda)
```

The long-range term keeps the Debye-Hückel physics — it must, since the limiting
law is exact — in the form

```math
f^{\gamma} = -A_\varphi\left[\frac{\sqrt{I}}{1 + b\sqrt{I}}
             + \frac{2}{b}\ln\left(1 + b\sqrt{I}\right)\right] ,
```

and the pair term carries the ionic-strength dependence that a single constant
cannot:

```math
B_{ca}(I) = \beta^{(0)}_{ca}
          + \beta^{(1)}_{ca}\,g\!\left(\alpha_1\sqrt{I}\right)
          + \beta^{(2)}_{ca}\,g\!\left(\alpha_2\sqrt{I}\right) ,
\qquad
g(x) = \frac{2\left[1 - (1+x)e^{-x}\right]}{x^{2}} .
```

``A_\varphi`` is the same Debye-Hückel slope as §1, on the osmotic basis:
``A_\varphi = A\ln 10/3``, so it comes from the water model rather than from a
second implementation. ``\alpha_1``, ``\alpha_2`` and ``b`` are not fitted —
they fix the functional form, and a published ``\beta`` table is only valid with
the values it was fitted against.

### Why the whole construction matters here

``\gamma_i`` and ``\varphi`` are both partial derivatives of **one** function.
So the Gibbs-Duhem relation between the solutes and the solvent is an identity
of the algebra, not an approximation — and that is the defect §3 measured in the
other two models. It is verified rather than asserted: with the derivative taken
analytically, ``\sum_i n_i\,\mathrm{d}\mu_i`` comes out **exactly zero** for
this model along three composition directions at 0.1, 1 and 3 mol/kg, where the
B-dot model's residual is not small.

### What it costs

The parameters are **caller input**, and there is no way around it: no
thermodynamic database ships them, and this package refuses to invent them. A
[`PitzerParameters`](@ref) takes every table as a keyword without a default, and
completeness is checked against the species list rather than against the set
alone — every cation-anion pair present must have a ``\beta^{(0)}``, and the
error names those that do not.

Two limitations belong here rather than in a footnote:

  - the **higher-order electrostatic terms** ``{}^E\theta(I)``,
    ``{}^E\theta'(I)`` are not implemented. They vanish identically for a
    symmetrical pair, so a single 1-1 or 2-2 electrolyte is unaffected; in a
    mixture of Na⁺ with Ca²⁺ — a cement pore solution — they are a real
    omission;
  - a published set is fitted **at one temperature**, and nothing here
    extrapolates the interaction parameters away from it.

### The speciation a set assumes is part of the set

This one is easy to miss and changes results. A Pitzer set absorbs ion
association into its ``\beta`` coefficients: Harvie, Møller and Weare fit
Ca–SO₄ interaction rather than postulating a `CaSO₄⁰` complex. A species list
that carries both the free ions **and** the ion pairs therefore counts the same
association twice, and CEMDATA18 does carry them — `Ca(SO4)@`, `CaOH+`,
`Na(SO4)-`, `NaOH@` among others. It also names the silica species differently
(`HSiO3-`, `SiO2@`) from the set's `H3SiO4-` and `H2SiO4-2`.

So the shipped Reardon set does not drop into a CEMDATA18 cement calculation,
and the completeness check refuses such a system rather than returning a number.
That refusal is the correct outcome and not a limitation of the implementation:
combining a dissociated parameterization with an associated speciation is not a
defensible calculation in any code.

## 7. Outside the domain

None of this survives arbitrary concentration, and the failure is silent: the
formulas go on returning finite, plausible numbers. Two guards exist rather than
one warning.

[`solvent_fraction`](@ref) reports the mole fraction of water *within* the
aqueous phase, and [`SOLVENT_FRACTION_FLOOR`](@ref) is where
[`equilibrate_certified`](@ref) refuses to call the answer a solution at all.
The regime it catches is real: a cement paste below ``w/c \approx 0.30``, on the
species lists used here, drives the free water to the solver's floor and the
ionic strength to hundreds of mol/kg — a composition at which every molality is
per kilogram of a solvent that is no longer there. See
[the w/c example](@ref sec-wc-ratio) and
[the water budget](@ref sec-theory-water-budget).

A Pitzer model does **not** rescue that regime. Its fitted range is a few
mol/kg, not hundreds, and no activity model repairs a composition that has left
the physical picture of a solution.
