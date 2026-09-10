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

```@example am
using ChemistryLab
using DynamicQuantities
using Printf

substances = build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false)
dict = Dict(symbol(s) => s for s in substances)
cs = ChemicalSystem([dict[s] for s in split("H2O@ H+ OH- Na+ Cl-")],
                    ["H2O@", "H+", "Na+", "Cl-", "Zz"])
sym_w = symbol(cs.species[only(cs.idx_solvent)])

function nacl(m)                       # m mol NaCl per kg of water, imposed
    st = ChemicalState(cs)
    set_quantity!(st, "H2O@", 1.0u"kg")
    set_quantity!(st, "Na+", m * u"mol")
    set_quantity!(st, "Cl-", m * u"mol")
    set_quantity!(st, "H+", 1.0e-7u"mol")
    set_quantity!(st, "OH-", 1.0e-7u"mol")
    return st
end
nothing # hide
```

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
[Helgeson1981](@cite) Table 1:

```@example am
for T in (298.15, 333.15, 373.15)
    p = hkf_debye_huckel_params(T, 1.0e5)
    w = water_thermo_props(T, 1.0e5)
    e = water_electro_props_jn(T, 1.0e5, w)
    ρ = w.D / 1000
    @printf("T = %6.2f K   ρ = %.4f g/cm³   ε = %6.2f   A = %.4f   B = %.4f\n",
            T, ρ, e.epsilon, p.A, p.B)
end
```

Both rise with temperature, because water's dielectric constant falls faster than
``T`` rises: hot water screens worse, so the same ionic strength costs more.

### The screening length is the size of a gel pore

The Debye length is ``\kappa^{-1} = 1/(B\sqrt{I})`` in ångström when ``B`` is in
Å⁻¹(kg/mol)^½, and it is worth putting a number on it:

```@example am
B25 = hkf_debye_huckel_params(298.15, 1.0e5).B
println("     I (mol/kg)    Debye length (nm)")
for I in (0.001, 0.01, 0.1, 0.3, 1.0, 3.0)
    @printf("   %10.3f    %14.3f\n", I, 1 / (B25 * sqrt(I)) / 10)
end
```

A cement pore solution sits around ``I \approx 0.1``–``0.5 mol/kg``, so its
screening length is a **few ångström** — the thickness of two or three water
molecules. That is the same scale as the water films in the gel pores of C-S-H,
and it is worth noticing that both assumptions of step 1 and step 2 above are
strained there: a continuum of bulk permittivity, and ions independent in a mean
field. Nothing in the formulas announces it. It is the physical reason to treat
an extended Debye-Hückel model as a correlation valid in bulk solution rather
than as a theory of confined water.

## 2. The ``\dot{B} I`` term is a deviation function, not a physical term

Measured activity coefficients turn back *upwards* at high ionic strength, which
no screening argument produces. The B-dot model adds a linear term for it,

```math
\log_{10}\gamma_i = -\frac{A z_i^2\sqrt{I}}{1 + B\mathring{a}_i\sqrt{I}}
                    + \dot{B} I .
```

Its status is worth being precise about. [AndersonCrerar1993](@cite) (§17.7.1,
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

| | [`DiluteSolutionModel`](@ref) | [`DaviesActivityModel`](@ref) | [`HKFActivityModel`](@ref) |
|:--|:--|:--|:--|
| solute scale | molarity | molality | molality |
| ``\gamma_i`` | ``\equiv 1`` | Davies | extended D-H + ``\dot{B} I`` |
| ``a_w`` | Raoult | Raoult | osmotic coefficient |
| per-species data | none | none | ion radii ``\mathring{a}_i`` (tabulated, overridable) |
| scalar inputs | none | ``A``, ``b``, ``b_n`` | ``A``, ``B``, ``\dot{B}``, ``K_n``, ``\mathring{a}_{\text{default}}`` |
| ``T``, ``P`` dependence | none | ``A(T,P)`` on request | ``A(T,P)``, ``B(T,P)`` on request |
| returns | ``\ln a_i`` for every species | same | same |
| ``\gamma`` useful to | ``I \lesssim 0.01`` | ``I \lesssim 0.5`` | ``I \lesssim 1`` |
| Gibbs-Duhem consistent | approximately | **no** (§5) | to ``10^{-5}`` (§5) |

All three return the same object — a vector of ``\ln a_i`` indexed like
`cs.species`, covering solutes, solvent, pure crystals (``0``), gases and
solid-solution end-members — so they are interchangeable at every call site, and
`concentration_scale` tells the accessors which convention was used.

## 5. Measured: where the three part company

Nothing here is solved. The models are evaluated on the same imposed NaCl
composition, which is the cheapest way to see what the choice is worth.

```@example am
models = ["dilute" => DiluteSolutionModel(),
          "Davies" => DaviesActivityModel(),
          "B-dot" => HKFActivityModel()]

println("               γ(Na⁺)                        a_w")
println("  m      dilute   Davies    B-dot      dilute   Davies    B-dot")
for m in (0.001, 0.01, 0.1, 0.5, 1.0, 3.0)
    st = nacl(m)
    γ = [activity_coefficients(st, mod)["Na+"] for (_, mod) in models]
    aw = [exp(log_activities(st, mod)[sym_w]) for (_, mod) in models]
    @printf("%6.3f  %7.4f  %7.4f  %7.4f    %7.5f  %7.5f  %7.5f\n", m, γ..., aw...)
end
```

Read the ``\gamma`` columns first. The ideal model is already several percent off
at a **millimolal**, which is worth knowing before treating ideality as a safe
default. The two corrections do not agree with each other either — they part
company around a tenth molal — and by 3 mol/kg Davies has returned
``\gamma > 1`` while the B-dot model is still below 1: the ``bI`` term has taken
over completely, which is the ceiling of §2 arriving.

Now the ``a_w`` columns, and here the surprise: they barely separate at all. The
Raoult and osmotic routes differ by a few parts in a thousand even at 3 mol/kg.
It would be easy to conclude that the water-activity route does not matter.

It does, and the next table is why.

### Gibbs-Duhem: the values agree, the derivatives do not

Equilibrium is set by chemical potentials, that is by *derivatives* of the
activities with respect to composition, not by their values. So the test that
matters is whether ``\sum_i n_i \,\mathrm{d}\mu_i = 0`` holds along a
composition change. It is measured here along three different directions,
because each exposes a different defect:

```@example am
using LinearAlgebra

M_W = 0.0180153
n_w = 1.0 / M_W
cs3 = ChemicalSystem([dict[s] for s in split("H2O@ Na+ Cl-")], ["H2O@", "Na+", "Cl-"])

function gd_residual(mod, m, dn)
    μ = build_potentials(cs3, mod)
    p = (ΔₐG⁰overT = zeros(3), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    n0 = [n_w, m, m]
    δ = 1.0e-6
    dμ = (μ(n0 + δ * dn, p) - μ(n0, p)) / δ
    return abs(sum(n0 .* dμ)) / max(norm(n0 .* abs.(dμ)), 1.0)
end

for (name, dn) in ("dissolution   dn = (0, +1, +1)" => [0.0, 1.0, 1.0],
                   "ion exchange  dn = (0, +1, -1)" => [0.0, 1.0, -1.0],
                   "water removal dn = (-1, 0, 0)" => [-1.0, 0.0, 0.0])
    println("\n── ", name)
    println("   m         dilute       Davies        B-dot")
    for m in (0.1, 0.3, 1.0, 3.0)
        r = [gd_residual(mod, m, dn) for (_, mod) in models]
        @printf("%6.2f    %10.3e   %10.3e   %10.3e\n", m, r...)
    end
end
```

Three readings, and they are the substance of this page.

**Along a true dissolution**, the B-dot model is four orders of magnitude more
consistent than Davies. And Davies is **worse than assuming ideality** — which
is not a paradox but the direct consequence of §3: correcting the solutes while
leaving the solvent at ``a_w = x_w`` makes the two halves of the model actively
contradict each other, whereas the ideal model at least contradicts itself less.
A model can be *more* wrong for being *partly* corrected.

**Along an ion exchange** at constant ``I`` and constant ``\sum m``, Davies and
the ideal model are indistinguishable — their coefficients depend on ``I``
alone, which does not move — and the only residual left is the B-dot model's own
defect: the charge-weighted mean radius of §3, showing up at a few parts in a
thousand. This is the direction `test/activities.jl` uses, which is why its
tolerance is `5e-3` and not solver tolerance.

**Along water removal** — the direction a drying paste actually takes — the
ordering is the same as for dissolution, and the gap widens as the solution
concentrates.

So the water-activity route is not a refinement on a number that hardly moves.
It decides whether the model is one thermodynamic system or two halves that
disagree, and that is visible only in the derivatives.

## 6. Outside the domain

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
[the w/c example](@ref sec-wc-ratio).

What is missing from this chapter is an **ion-interaction model** — Pitzer-class
— which replaces the correlation of §2 by a virial expansion of the excess
Gibbs energy: one coefficient per ion pair, one per triplet, and ``\gamma`` and
``\varphi`` derived from the same function so that Gibbs-Duhem holds by
construction. [AndersonCrerar1993](@cite) (§17.8) derive it as a cluster
expansion with osmotic pressure in place of pressure, which is what makes it a
different kind of object from anything on this page.
