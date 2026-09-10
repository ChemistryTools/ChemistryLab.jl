# [Thermochemistry: the quantities and how they connect](@id sec-theory-thermo)

This page is the shortest path from "what is a chemical potential" to the
expressions this package actually evaluates, with the **notation of the code**
rather than of a textbook. Everything here is either a definition or an identity;
nothing is fitted.

## 1. One potential, and the two pieces it is made of

The Gibbs energy of a system holding `nᵢ` moles of each species is
``G = \sum_i n_i \mu_i``, and the chemical potential of a species is what one
more mole of it costs:

```math
\mu_i = \left(\frac{\partial G}{\partial n_i}\right)_{T,P,n_{j\neq i}}
      = \mu_i^\circ + RT\ln a_i .
```

Two pieces, and they come from two different places in the package:

  - ``\mu_i^\circ(T,P)``, the **standard potential**, is thermodynamic data — it
    belongs to the species alone and comes from the database (§3);
  - ``a_i``, the **activity**, is where the mixture enters, and it comes from the
    [activity model](@ref sec-theory-activity).

Divide by ``RT`` and the code appears verbatim. [`build_potentials`](@ref)
returns the closure

```julia
μ(n, p) = p.ΔₐG⁰overRT .+ lna(n, p)      #  μᵢ/RT  =  ΔₐG⁰ᵢ/RT  +  ln aᵢ
```

so throughout this documentation

```math
g_i \;\equiv\; \frac{\mu_i}{RT} \;=\; \underbrace{\texttt{ΔₐG⁰overRT[i]}}_{\text{database}}
      \;+\; \underbrace{\texttt{lna(n,p)[i]}}_{\text{activity model}} .
```

Working in ``\mu/RT`` is not cosmetic: it is dimensionless, it is what the
minimizer differentiates, and it makes ``\ln a_i`` and the standard potential
directly comparable in size.

## 2. Activities, standard states, and the conventions in use

An activity is always a ratio to a **standard state**, and the choice of
standard state is a convention that must be stated or the numbers are
meaningless. The three in use here:

| class | activity | standard state | convention |
|:--|:--|:--|:--|
| aqueous solute | ``a_i = \gamma_i\, m_i/m^\circ`` | ``m^\circ = 1`` mol/kg of solvent | molality |
| aqueous solute (ideal model) | ``a_i = c_i/c^\circ`` | ``c^\circ = 1`` mol/L | molarity |
| solvent (water) | ``a_w`` | pure water | mole fraction |
| pure solid, pure phase | ``a_i = 1`` ⇒ ``\ln a_i = 0`` | the pure substance | — |
| gas in an ideal mixture | ``a_i = x_i`` | the pure gas at ``P`` | mole fraction |
| solid-solution end-member | ``a_k = \gamma_k x_k`` | the pure end-member | mole fraction (§5) |

[`concentration_scale`](@ref) is how a model declares which of the first two it
uses, and the aqueous accessors read it rather than guessing:

```@example thermo
using ChemistryLab
concentration_scale(DiluteSolutionModel()), concentration_scale(HKFActivityModel())
```

An activity coefficient is therefore not a property of a species but of a
species **in a solution**, and ``\gamma_i \to 1`` in the standard state's own
limit — infinite dilution for a solute, purity for a solid.

## 3. Where ``\mu^\circ(T,P)`` comes from

The database stores, per species, the coefficients of a heat-capacity
polynomial, and the package integrates it. That is the whole of the temperature
dependence, and it is exact rather than a series truncation: `THERMO_MODELS`
(`src/thermodynamics/thermo_models.jl`) holds ``C_p(T)`` and the *analytic
integrals* of the same coefficients for ``S``, ``H`` and ``G``, so the three
identities

```math
C_p = \left(\frac{\partial H}{\partial T}\right)_P ,
\qquad
S(T) = S(T_r) + \int_{T_r}^{T}\frac{C_p}{T}\,\mathrm{d}T ,
\qquad
G = H - TS
```

hold by construction rather than numerically. Aqueous solutes use a second
model, `:solute_hkf88_reaktoro`, the Helgeson-Kirkham-Flowers equation of state,
whose electrostatic part is the Born term evaluated with the dielectric constant
of water — the same water model that supplies ``A`` and ``B`` to the activity
models.

Two consequences worth carrying:

  - a reaction's temperature dependence is not an input; it follows from the
    ``C_p`` of its reactants and products, which is the **Gibbs-Helmholtz**
    relation ``\partial(\Delta_r G^\circ/T)/\partial T = -\Delta_r H^\circ/T^2``
    (van 't Hoff, when ``\Delta_r H^\circ`` is taken constant);
  - a species whose database entry has no ``C_p`` model cannot be moved off
    ``T_r`` at all, and the package says so rather than extrapolating.

## 4. Equilibrium is a constrained minimization, and its dual is the useful part

At fixed ``T`` and ``P``, equilibrium is the composition that minimizes ``G``
subject to the matter available. "The matter available" is a linear constraint,
because atoms and charge are conserved whatever the reactions do:

```math
\min_{n \ge 0} \; \sum_i n_i\, g_i(n)
\qquad\text{subject to}\qquad
\mathbf{A}\,n = b ,
```

where ``\mathbf{A}`` is the conservation matrix `cs.SM.A` — one row per
component, one column per species — and ``b`` the element (and charge) budget.
This is the formulation of [Leal2017](@cite), and it is why no reaction list is
needed: reactions are whatever moves ``n`` inside the null space of
``\mathbf{A}``.

The Lagrangian multipliers of the equality constraints are the interesting
output. Writing ``y_c`` for the multiplier of row ``c``, the first-order
conditions are

```math
g_s \;=\; \sum_c A_{cs}\, y_c \quad\text{for every species actually present},
\qquad
g_s \;\ge\; \sum_c A_{cs}\, y_c \quad\text{for every one absent.}
```

The ``y_c`` are **component potentials** (element potentials, when the components
are elements): one number per component, from which the chemical potential of
any species is a dot product. The two lines above are the complementarity
conditions of the KKT system, and they are exactly what
[`optimality_certificate`](@ref) checks — see
[Proving that an answer is the answer](@ref sec-theory-certificate) for the certificate itself and why
convexity makes it sufficient rather than merely necessary.

## 5. Reaction quotients, ``K``, and the saturation index

Take a species ``s`` and write its formation from the system's primary species,
``s = \sum_c A_{cs}\,(\text{primary }c)`` — the column of ``\mathbf{A}`` *is*
that reaction. For any reaction,

```math
\Delta_r G = \Delta_r G^\circ + RT\ln Q ,
\qquad
\Delta_r G^\circ = -RT\ln K ,
```

so ``\Delta_r G = RT\ln(Q/K)``: the sign of ``\ln(Q/K)`` says which way the
reaction runs, and ``Q = K`` is equilibrium. For a solid dissolving into its
ions, ``Q`` is the **ion activity product** IAP and ``K`` is the solubility
product, and the logarithm of their ratio is the **saturation index**.

Because the rows of ``\mathbf{A}`` are labeled by primary species, a component
potential *is* that primary's ``g``, and the index becomes a difference of
potentials with no equilibrium constant to look up:

```math
\mathrm{LogSI}_s \;=\; \frac{1}{\ln 10}\left(\sum_c A_{cs}\, y_c \;-\; g_s\right)
\;=\; \log_{10}\frac{\mathrm{IAP}}{K_{sp}}
\;=\; -\frac{\Delta_r G}{RT\ln 10} .
```

Negative means undersaturated, zero means in equilibrium with the solution,
positive means the phase **should have precipitated**. This is what
[`saturation_indices`](@ref) returns.

Two consequences of the identity are worth carrying, because they are what make
the index trustworthy rather than merely conventional:

  - **no equilibrium constant is looked up.** ``K_{sp}`` never appears in the
    computation; it is implied by the standard potentials, so an index cannot
    disagree with the ``\Delta_a G^\circ`` the rest of the calculation used;
  - **every phase present at an equilibrium must come out at exactly zero.**
    That is not a property of the answer, it is the first-order condition of
    §4 restated, so it is a check on the solve: measured on a CEM I paste, the
    twelve present solids land within ``1.2\times10^{-12}``. If they do not,
    the state is not an equilibrium and no other index in the result means
    anything.

The identity itself — that the difference of component potentials above equals
what [`saturation_indices`](@ref) computes — is asserted in the test suite
(`test/aqueous_properties.jl`), where it is recomputed by hand from this formula
and compared to the function.

## 6. Volume, porosity and chemical shrinkage

Volumes are treated as **ideal**: the volume of a phase is the sum of its
species' standard molar volumes,

```math
V = \sum_i n_i V_i^\circ(T,P) ,
```

with no excess volume of mixing. That is an assumption, and it is the one behind
every porosity this package reports. [`volume`](@ref) returns the split by
aggregate state, [`porosity`](@ref) the void fraction relative to a reference
state, and [`chemical_shrinkage`](@ref) the volume the reaction itself consumes —
hydrates occupying less than the water and clinker they were made from, which is
the mechanism behind self-desiccation
([Self-desiccation](@ref sec-self-desiccation)).

A species carrying no ``V^\circ`` contributes zero and would corrupt a porosity
in silence, which is why the volume machinery reports
`missing_molar_volumes(state)` and why [`CapillaryWater`](@ref) refuses to be
constructed when that list is not empty.

## 7. What is *not* assumed, and what is

The two lists, side by side:

| assumed | not assumed |
|:--|:--|
| one well-mixed phase per aggregate state | a reaction list, a reaction path, a sequence |
| ideal molar volumes, no excess volume | ideal activities — that is the model's business |
| the activity model's own domain of validity | that the answer is a local minimum: the certificate proves global (see below) |
| that ``\mathbf{A}n = b`` is the whole of the conservation | which phases appear: they are a result |

Convexity is what makes the last cell true. ``G`` is convex in ``n`` for the
built-in activity models, so a KKT point is *the* minimum and the certificate is
sufficient rather than necessary — the one documented exception being
[`CapillaryWater`](@ref), where a composition-dependent shift of ``\ln a_w``
breaks the convexity argument and the certificate falls back to proving a KKT
point.

## See also

  - [Activity models](@ref sec-theory-activity) — the ``\ln a_i`` half of ``g_i``
  - [Solid solutions](@ref sec-theory-solid-solutions) — the mole-fraction half
  - [Chemical Equilibrium](@ref sec-equilibrium) — driving the solver
  - [Proving that an answer is the answer](@ref sec-theory-certificate) — the certificate
