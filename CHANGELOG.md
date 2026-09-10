# Changelog

## v0.17.0 — the water that is there, and the water that counts

Two halves. An **ion-interaction activity model**, which is the second of the
two ingredients the roadmap named for predicting where a sealed paste stops; and
the **coupled kinetic run** that 0.16.0 described and never performed, which
turns out to have been impossible for a reason worth recording. Around them, the
documentation was reorganized into four chapters, and every activity model in
the package now states the physics behind its formulas and the provenance of
every default it carries.

### Breaking changes

- **`ΔₐG⁰overT` is renamed `ΔₐG⁰overRT`.** It was a typo: the quantity is
  `ΔₐG⁰/(RT)`, dimensionless, and it has to be, since it is added directly to
  `ln aᵢ` to form `μᵢ/RT`. The parameter tuple is part of the documented
  interface — an `activity_model` closure reads `p.ΔₐG⁰overRT` — so a custom
  activity model written against the old name now fails with a
  `NamedTuple has no field` error rather than silently. 45 occurrences renamed.
- **The registry treats a minor bump below 1.0 as breaking whatever the API
  did**, so `[compat] ChemistryLab = "0.16"` will not accept `0.17`. Downstream
  packages must widen their bound. `MeanFieldHomogenization.jl` depends on this
  package only in `docs/Project.toml`.
- **`PoreHumidity` no longer differentiates its retention law at full
  saturation.** Its value there is unchanged; its derivative is now zero instead
  of infinite. Anything that read the derivative at `S = 1` — nothing could,
  since the integration it exists for did not run — changes.

### Added — `PitzerActivityModel`

An ion-interaction model, and a different kind of object from everything else in
`activities.jl`. Those are corrected Debye-Hückel laws: one screening term, one
size correction, one empirical term for the rest. This one expands the excess
Gibbs energy as a **virial series in the molalities** — a coefficient per ion
pair, one per triplet, no per-species radius — which Anderson & Crerar derive as
a cluster expansion with osmotic pressure in place of pressure.

The consequence is the reason to have it: `γ` and the osmotic coefficient are
partial derivatives of **one** function, so the Gibbs-Duhem relation between
solutes and solvent is an identity of the algebra. Measured, with the derivative
taken analytically: the residual is **exactly zero** along three composition
directions at 0.1, 1 and 3 mol/kg, where the B-dot model's is not small. A
finite difference cannot see this — its own truncation error at 0.1 mol/kg is
larger than the quantity — which is why the test uses AD.

Also validated: the dilute limit reproduces `log₁₀ γ± → −A|z₊z₋|√I`, and the
model agrees with the independent B-dot implementation at a millimolal, as two
models sharing a limit must.

**Nothing is parameterized by default, and completeness is checked against the
species list rather than against the set.** A table is complete or not
*relative to a system*, so `PitzerParameters` takes every table as a keyword
without a default — `UndefKeywordError` before a number is computed — and the
check that matters happens when the model meets a `ChemicalSystem`: every
cation-anion pair present must have a `β⁰`, and the error names those that do
not. A missing `θ`, `ψ` or `λ` is zero, which is the convention of the
literature the tables come from; a missing `β⁰` cannot be, because falling back
on ideal behavior for one pair of a Pitzer calculation is not an answer.

Two limitations are in the docstring rather than left to be discovered: the
higher-order electrostatic terms are not implemented, which is exact for a
symmetrical pair and an omission in a Na/Ca mixture; and no temperature
dependence of the interaction parameters themselves.

### Added — a cited cement parameter set, with its provenance per entry

`data/pitzer-reardon1990.toml`, read through `build_pitzer_parameters`, from
Reardon (1990), *Cement and Concrete Research* **20**, 175–192, whose tables are
those of Harvie, Møller & Weare (1984) with three exceptions he names. Nothing
about it is automatic: the file has to be named.

Two properties of the set change what a user should conclude from it, so they
are recorded in the file and readable from code:

- the silicate, aluminate and ferrate parameters are **estimates, not
  measurements**. Reardon says which analog each borrows — HSO₄⁻ for Fe(OH)₄⁻,
  Al(OH)₄⁻ and H₃SiO₄⁻, SO₄²⁻ for H₂SiO₄²⁻, H₂CO₃⁰ for H₄SiO₄⁰ — and those are
  exactly the ions a cement assemblage needs. Every affected entry carries
  `origin = "estimated:<analog>"`, read back by `pitzer_origin`;
- the set assumes a **fully dissociated** speciation. The association of Ca with
  SO₄, of Na with OH, is already inside these β coefficients, so a species list
  that also carries `Ca(SO4)@` or `CaOH+` — as CEMDATA18 does — counts each
  association twice, and CEMDATA18 also names the silica species differently.
  The completeness check refuses such a system rather than returning a number,
  and that refusal is the correct outcome in any code, not a shortcoming here.

Transcription was verified rather than trusted: Table 2's 216 values were
checked against the PDF text layer and all matched, the only surplus tokens on
the PDF side being OCR fragments of the Greek symbols. Tables 5 and 6 could not
be checked that way — that text layer turns one −0.0677 into −0.0077 and drops
two others — so they were read from the pages rendered at six times
magnification, and the data file's header says so.

### Fixed — `PoreHumidity` could not be integrated

0.16.0 built it, wired the dispatch, and said in the self-desiccation page that
handing it to `parrot_killoh_avrami` and integrating was how to obtain the
arrest in time. Nothing ever did, and doing it found why: **the integration did
not advance past `t = 0`.**

Every retention law of the van Genuchten family has an unbounded `dp_c/dS` at
full saturation, and a sealed paste starts saturated. An implicit solver
differentiates the rate law at its first step, so it received an infinite
Jacobian entry — the gradient of the humidity with respect to the composition
came back with an infinite norm — and returned immediately with every degree of
hydration at zero. No unit test on the object could see this: the object is well
behaved, and only an integration reaches the singular point.

The value at saturation is not in doubt, so it is now returned as a constant and
the derivative there is zero. The code says plainly that this is a
**regularization and not an identity**, that it applies within `1e-10` of
saturation — the initial condition alone, in practice — and that the solver
leaves that point on its first successful step, after which the real and stiff
derivative applies.

What the coupled run then gives, with `α_max = 1.0` and no Powers bound: the
arrest is a **result**. The uncoupled rate reaches the same α(C3S) = 0.9164 at
every w/c, because it knows nothing about how much water there is; coupled, α
runs from 0.556 at w/c 0.25 to 0.837 at 0.50. The drier mixes stop at an
internal humidity of 0.800 and a pore saturation of **0.786** — which is the
`S*` the static water budget reads off the same measured retention curve at
RH 0.80. The trajectory arrives at the number the budget assumes, by a different
route, with nothing arranged to make it so.

It also shows that the proportionality `α_max ∝ w/c`, exact in the static
construction, does **not** survive integration: `k = (w/c)/α_max` is near 0.45
at the dry end and rises with w/c. The static page is careful to say that the
proportionality is structural and no evidence; this is what it looks like when
the structure is removed.

### Added — `water_activity`, and a `kelvin_shift` on `log_activities`

0.16.0's own release notes recorded that `log_activities` did not know about the
capillary shift, so a solve posed at `a_w = 0.90` read back as 0.999995. Both
accessors now take a `kelvin_shift` keyword, `0.0` by default so that nothing
changes for a caller who does not ask, and `water_activity(state, model;
kelvin_shift)` returns the composed value. The chemical and capillary lowerings
are independent, their chemical potentials add, and so their activities
multiply — which is what the keyword implements and what the docstring derives.

### Documentation — four chapters, and the physics behind every formula

The page tree was four flat lists in the order the pages were written. It is now
grouped, and each chapter answers one question: **Theory** why this is the right
calculation, **Manual** how an object is written, **Tutorials** how to drive a
calculation, **Applications** what a real case looks like and what the choices
cost in numbers. Nine pages that were object syntax filed under Tutorials moved
to the Manual. **Cementitious media is a subsection of four chapters rather than
a chapter of its own**, because the material is the subject of the package and
belongs wherever its question is being asked. No page was removed; three paths
that released notes link by URL were deliberately left where they are.

A **Theory** chapter, with the rule written on its own index that a theory page
may show code where the correspondence with the implementation is the point, and
does not build systems, solve, or print tables:

- **Thermochemistry** — the definitions and identities in the code's own
  notation, where `μ°(T,P)` comes from, equilibrium as a constrained
  minimization with the component potentials as its multipliers, and the
  saturation index as a difference of those potentials. That identity is now
  asserted in the test suite, which recomputes `LogSI` by hand from the formula
  the page states.
- **Proving that an answer is the answer** — the convexity proof, the KKT
  conditions and the dual solver, moved out of the tutorial they were living in.
- **Activity models** — where the `√I` comes from, in three steps with the
  assumptions each makes, since those assumptions are what later fails; why `A`
  and `B` are properties of water rather than fitting constants; the Debye
  length, which at a cement pore solution's ionic strength is 0.55 nm — the
  width of the gel pores where a paste keeps its last water, so the
  continuum-dielectric and mean-field assumptions are strained exactly where the
  arrest happens; and `Ḃ` in its real status as a deviation function fitted to
  one salt.
- **Solid solutions** — mixing entropy, the excess Gibbs energy, the Margules
  and Redlich-Kister forms as implemented, and what the sign of `W` means. Two
  identities are checked rather than claimed: that Redlich-Kister reduces to a
  regular solution, and that `W = 2RT` is the critical point.
- **The water budget of a hydrating paste** — why a Gibbs minimization predicts
  that clinker survives below `w/c ≈ 0.30` while Powers reports 0.42, and why
  the gap is not a thermodynamic quantity: the capillary term is two orders of
  magnitude too weak, and what stops a real paste is a broken liquid path
  across four orders of magnitude of length, which a 0D model has no
  representation of.

Every default of every activity model is now classified as derived, tabulated,
or **a convention with no source recorded in this package** — `Ḃ = 0.041`,
`Kₙ = 0.1` and `å_default = 3.72` are in the third class, and that is said out
loud rather than given a plausible citation.

Long solver output is folded into collapsed blocks rather than unrolled: seven
pages ended on an assignment whose value Documenter then displayed in full, up
to 48 species and a whole conservation matrix between two paragraphs.

### Compatibility

Tested on Julia 1.12 and 1.13. `[compat] julia = "1.12"` is unchanged and
already admitted 1.13.


## v0.16.0 — the water a paste cannot use

A sealed cement paste stops hydrating before it runs out of cement, and until now
nothing in the package represented why. `powers_alpha_max` supplied the empirical
cap and said in its docstring that the bound is "a statement about transport and
access, not about thermodynamics". This release supplies the missing state
variable — the internal relative humidity of a partially saturated pore space —
and the tutorial that validates it recovers Powers' coefficient from independent
data.

### Breaking changes

Nothing in the API breaks and no default behavior changes. Two consequences are
breaking in practice:

- **The registry treats a minor bump below 1.0 as breaking whatever the API did**,
  so `[compat] ChemistryLab = "0.15"` will not accept `0.16`. Downstream packages
  must widen. `MeanFieldHomogenization.jl` depends on this package only in
  `docs/Project.toml` and needs `"0.14, 0.15, 0.16"`.
- **`solve_certified` now writes the *winning* candidate's parameters** into the
  caller's `parameters` Ref. It wrote the last candidate's while returning the
  best-by-error one, so a prescribed-pH scan reported a titrant amount belonging
  to a different composition. Code that relied on the old value was reading a
  mismatch.

### Added — water retention, and the Kelvin relation

`WaterRetention` and its subtypes describe how tightly a pore space holds the
water still in it. That is a constitutive input, measured on a particular material
at a particular age, so **nothing ships with a default value**: every named law
takes its parameters as keywords without defaults, and Julia raises
`UndefKeywordError` before any number is computed. The enforcement is the
language's, not a check that can be forgotten.

- `TabulatedRetention(; S, a_w)` takes a measured isotherm and interpolates
  linearly in `ln a_w`, the variable the coupling reads. Its inner constructor
  refuses a non-monotone table rather than interpolating nonsense, and leaving the
  measured range clamps and warns once instead of extrapolating a logarithm into
  confident absurdity.
- `VanGenuchten(; a, m)`, with the published parameters of Baroghel-Bouny et al.
  (1999) quoted in its docstring together with the convention they use — they
  write the same expression with `b = 1/m`, and getting that inversion wrong is
  silent.
- `FunctionRetention(f)`, the escape hatch: a bare function of the saturation,
  returning the activity directly. It is the one law nothing validates, since
  `f` is opaque, and its docstring says so — `CapillaryWater` still tests it at
  `S = 1`. `CapillaryWater` and `PoreHumidity` wrap a bare function in it
  themselves, so no caller has to name it.
- `kelvin_activity`, `kelvin_radius`, `capillary_pressure`, `water_activity`.
  RH 80 % is a meniscus of radius 4.8 nm, the gel-pore scale, which is why the
  water Powers assigns to gel pores and the water a sealed paste cannot use are
  the same water.

### Added — `CapillaryWater`, a constraint with a certificate

One unknown, the Kelvin shift of the solvent's chemical potential; one equation,
the retention law at the current saturation `S = V_liquid/(V_ref − V_solid)`. No
titrant column, because a sealed specimen exchanges water with nothing — the whole
difference from `FixedActivity`. It refuses eagerly when a species present carries
no standard molar volume, since such a species contributes zero to the pore volume
in silence and would corrupt the saturation with nothing to show for it.

Measured on a CEM I paste: certified at every w/c from 0.30 to 0.50, constraint
residual between 1e-16 and 1e-13, internal humidity falling from 0.87 to 0.32 as
the mix gets drier.

Two caveats are in the docstring rather than left to be discovered. The
certificate proves a **KKT point of the constrained problem** — stationarity, mass
balance, no absent phase supersaturated, and the capillary closure satisfied — but
not global optimality, because a composition-dependent shift of `ln a_w` is not
derived from a convex `G` for an arbitrary retention law. And **`log_activities`
does not know about the shift**: it evaluates the activity model, so with the
shift at `ln 0.90` it returns `a_w = 0.999995` while the solve was posed at 0.90.

### Added — `PoreHumidity`, and where the arrest actually comes from

The Kelvin term is **two orders of magnitude too weak to arrest hydration**.
Measured, with a certificate on every answer: imposing a water activity from
saturation down to 0.80 leaves the equilibrium assemblage of a CEM I paste
unchanged to six digits. Below that the constrained problem stops certifying, so
nothing is claimed there — and the arithmetic says nothing is expected there
either. At `a_w = 0.80` the shift is
553 J per mole of water, worth about 1.8 kJ per mole of alite against a hydration
Gibbs energy of order −100 kJ/mol; nulling it would need `a_w ≈ 5e-6`, a Kelvin
radius smaller than a water molecule. A real paste stops at 75–80 % RH because
transport and nucleation stop.

So the humidity belongs in the rate law, and `humidity_factor` has implemented
that cut all along. What was missing is that its argument was exogenous:
`_humidity_at` took a constant or a function of **time**. The rate closure already
receives the full composition, so `PoreHumidity(retention, system; reference)`
computes the saturation and returns the pore humidity, and a third method reads
it. The other two ignore the new argument, so nothing a caller wrote before
changes.

### Fixed — the certificate now audits the problem that was solved

`optimality_certificate` rebuilt a `FixedTP` problem whatever constraint had been
applied. That misses two things at once: the conservation rows lose their `Aq q`
term, so a prescribed activity or pH is measured against a budget short by exactly
the titrant amount; and `hq` is skipped, so a constraint that shifts a chemical
potential is measured against the unshifted one and can never certify however
right it is. It takes `constraint` and `q` keywords, both defaulting to the old
behavior, and returns `param_residual`.

`_kkt_error` counts that residual, so the multi-start search can no longer prefer
a candidate that minimizes the Gibbs energy while violating the very equation that
makes it a constrained answer.

### Compatibility

Tested on Julia 1.13, which the CI matrix now runs explicitly alongside the
1.12 floor. `[compat] julia = "1.12"` already admitted it and did not change.

### Documentation — Powers' 0.42, taken apart

A new tutorial, [Self-desiccation: where Powers' 0.42 comes from](https://micropochemomechanics.github.io/ChemistryLab.jl/stable/tutorials/self_desiccation/),
with `scripts/self_desiccation_powers.jl` as its runnable companion.

It derives `α_max = (w/c)/k` with `k = b + s S*/(1−S*)`, then fills each term from
its own source: `b = 0.3095` g/g and `s = 0.0639` cm³/g from certified equilibria
at imposed degree of hydration (both constant in α to four digits, which is the
check on the linear budget the derivation assumes), and `S* = 0.7862` at RH 0.80
from a measured desorption isotherm.

**Inverted, Powers' 0.42 corresponds to an arrest at 77.5 % relative humidity** —
assembled from a chemical shrinkage computed out of the thermodynamic database, an
isotherm fitted for a drying study two decades earlier, and Powers' own water
split. That is the window sealed pastes are independently reported to stop in.

The page is explicit that the *proportionality* `α_max ∝ w/c` is structural and
therefore no evidence, that only the coefficient is predicted, and that the
residual 11 % is attributable: the model's formula water exceeds Powers'
non-evaporable water by 0.0795 g/g, which is the interlayer water CEMDATA18 writes
into the C-S-H formula and D-drying removes. Not the same quantity, so not an
error — and the `CSHQ` solid solution binds more, 0.3684 g/g, so it moves away
from Powers rather than toward him.

It also **executes the two controls that could have refuted it**, rather than
asserting them: the imposed-activity sweep above, which is why the arrest is read
from the rate law and not from the Gibbs energy; and the same construction run on
a second measured curve from the same table, which returns a ratio
`α_max/(w/c)` constant to every digit printed in both cases and differing between
them by 18 %. The constancy is an identity of the construction and carries no
information; the value, `1/k`, carries all of it.


## v0.15.2 — a solution that is no longer a solution

A certificate proves that a composition minimizes the Gibbs energy of the problem
**as posed**. It says nothing about whether the problem was posed inside the
model's domain, and there is one way to leave that domain while producing an
ordinary-looking `ChemicalState`: let the solids take all the water.

### Added — `solvent_fraction`, and a guard on the answer

`solvent_fraction(state)` is the mole fraction of the aqueous solvent within the
aqueous phase, `n_w / Σ_aqueous n`. Every quantity built on that phase —
molality, ionic strength, activity, pH — is defined *per kilogram of solvent*,
and `DualEquilibriumSolver` parameterizes its interior variables by the solvent's
chemical potential. All of it presumes the solvent **is** the phase. A real
electrolyte keeps `x_w` above about 0.9: seawater is 0.99, a saturated NaCl brine
0.90.

Measured on a sealed cement paste taken below its stoichiometric water demand —
w/c = 0.28 on the mix of the w/c example — the Gibbs minimum consumes the free
water down to **6e-9 mol**, the solver's floor. The solvent then holds a **fifth**
of its own aqueous phase and the ionic strength is reported as **409 mol/kg** by
a Debye-Huckel model valid to about one. Nothing in the certificate objects,
because nothing is wrong with the minimization: forming more hydrate always
lowers the energy, and nothing in the model penalizes a solution concentrated
past any physical meaning. The water activity of a real paste collapses as the
pores empty and stops the reaction; an activity model extrapolated that far goes
on returning finite numbers.

`equilibrate_certified` now checks `solvent_fraction` on the answer it returns
and says so, naming what happened, when it falls below `SOLVENT_FRACTION_FLOOR`.
Warned rather than raised under the default flag, because the composition of the
*solids* still carries the mass-balance information a caller may legitimately
want; under `STRICT_CONVERGENCE[]` it raises, like any other answer that is not
one.

`SOLVENT_FRACTION_FLOOR = 0.5` is not a modeling choice but a floor no real
solution comes near — about 28 mol of solute per kilogram of water, past
saturation for anything. A value below it means the solve drove the water into
the solids, not that the solution is concentrated.

### Documentation — a usable answer below the stoichiometric demand

A high-performance concrete is mixed at w/c between 0.25 and 0.35, which is an
ordinary regime and must be computable, with a certificate and with all three of
the things one wants from it: how much clinker stays unhydrated, which hydrates
form, and what the pore solution contains. The w/c example now works that case.

The construction is to stop the reaction where the physics stops it rather than
ask the minimizer to discover an arrest point it has no term for: react a fraction
α of the clinker with **all** the water, and leave the rest unhydrated. The
equilibrium is then computed on a system that still has a solution in it, and α
is what `powers_alpha_max` supplies. Imposing the reacted fraction is the standard
construction of cement thermodynamic modeling — it is how Lothenbach & Winnefeld
(2006) compute a hydrating paste.

Measured: at w/c = 0.25 with α = 0.595 the solve certifies, the solvent holds
0.9993 of its phase, the ionic strength is 0.0351 mol/kg and the pH 12.39, with
40.5 % of the clinker unhydrated and a total porosity of 0.2056 — against the
409 mol/kg and vanished solvent of the unconstrained solve at the same w/c.

Predicting the arrest point instead of imposing it is a well-posed thermodynamic
question the package cannot answer yet, and the page says which two ingredients
are missing. An **activity model valid at very high concentration**, Pitzer-class,
because what physically stops hydration is the collapse of the water activity as
the last of the pore solution is consumed, and an extended Debye-Huckel model
extrapolated to 409 mol/kg goes on returning finite numbers instead of collapsing.
And a **coupling between pore structure and water activity**, the Kelvin term,
because water in a fine pore is held at a reduced activity whatever its
composition — which is what self-desiccation is, and it is poromechanics rather
than solution chemistry.

### Documentation — the water-limited regime, and Powers' 0.42

`docs/src/examples/cement_wc_ratio.md` claimed that "no clinker survives at any
w/c", explained by the minimum "always forming a less hydrous assemblage rather
than leave alite standing". The first half is true only over the 0.30-0.60 range
it scans, and the second is false: the least hydrous assemblage available still
binds water, and when there is not enough, alite stays.

Measured on the page's own species list, clinker left as a fraction of the
cement: 54.6 % at w/c = 0.15, 34.3 % at 0.20, 14.1 % at 0.25, 2.0 % at 0.28, and
zero from 0.30 up. The crossover is a property of the hydrate assemblage, not a
constant, so the page measures it in an executed block instead of quoting one —
and reads it for one thing only, the *existence* of the regime, since past it the
aqueous phase has vanished and the amounts are not equilibrium values.

That does not put the page at odds with Powers, and the reason is worth stating
because the two numbers measure different things. Powers' 0.42 g of water per
gram of cement is **not** a stoichiometric demand: it is about 0.23 g of
non-evaporable water, written into the hydrate formulae and a mass balance any
Gibbs minimization must respect, plus about 0.19 g of **gel water** held in the
C-S-H gel pores, physically present and chemically unavailable. A sealed paste
stops by self-desiccation with water still in the specimen. So between roughly
0.23 and 0.42 the equilibrium and Powers disagree and **both are right**, and
below the stoichiometric demand they agree. `powers_alpha_max`'s docstring said
"hydrating one gram of cement binds about 0.42 g of water", which is the
misreading that produced the wrong claim; it now separates the two halves and
notes that curing water moves the bound to about 0.36.

The page also solves with `equilibrate_certified` rather than through Ipopt, and
prints whether every point certified. Ipopt is kept as a documented alternative
in an admonition, not executed: not a worse optimizer, but a bare interior point
that returns an iterate and no statement about it, and an extra binary dependency
for a calculation the package can already prove. The visible difference is small
and telling — absent phases come back at exactly zero instead of sitting at the
solver's 1e-8 lower bound. Every number on the page was re-measured on the
certified route and is unchanged: pH 12.3924 flat across the scan, porosity
12.2 % to 41.0 %, and at w/c = 0.50 the one-argument porosity 0.2846 against a
two-argument total of 0.3376 with the volume shrinking 7.41 %.


## v0.15.1 — a starting point that conserves mass, and a strict flag that is read

The automatic initial approximation shipped in 0.15.0 could hand the solver a
starting point that is not on the constraint surface, and then build every later
step on it. Reported from a Windows run of the same cement paste that certifies
here: `optimal = false`, an element balance off by **6.7 mol**, every hydrate at
zero — and a table of amounts that reads like a result.

### Added — the certificate names the missing phase, and the route puts it in

`saturation_indices(state, model)` returns `LogSI` for every species: zero for a
phase at equilibrium with the solution, negative for an undersaturated one,
positive for one that **should have precipitated**. It is what GEM-Selektor
prints as `LogSI`, and what `optimality_certificate` had been compressing into a
single worst violation without saying which phase it was.

No fitting is involved: the row labels of the conservation matrix are the primary
species, so a component's element potential is that primary's chemical potential
and `LogSI_s = [Σ_c A_cs μ_c/RT − μ_s/RT] / ln 10`. The check comes with it —
every phase actually present at an equilibrium must come out at zero, and on a
CEM I paste the twelve present solids land within 1.2e-12.

`equilibrate_certified` now **acts** on that. When its answer is uncertified
because a phase sits at the lower bound while the solution is supersaturated with
respect to it, the route puts that phase in and solves again, at most
`_MAX_RESTARTS` times.

The failure this exists for is a **phase swap**, which an active-set loop that
admits one phase at a time cannot perform. Reported from one machine while
another certified the same source: all 0.02515 mol of magnesium sat in brucite
with `hydrotalcite` absent and supersaturated by 5.58 log units, and admitting
the hydrotalcite requires dissolving the brucite entirely and taking aluminum
back from the hydrogarnet in the same step. The solve was otherwise impeccable —
stationarity 1.5e-16, element balance 1.8e-14 — which is exactly what an answer
converged onto the wrong active set looks like.

It is reproduced exactly, and now fixed: solving the paste with `hydrotalcite`
out of the phase list gives a certified 74.1514 cm3 with all the magnesium in
brucite, which is the reported state to seven digits; putting hydrotalcite back
and starting from there,`equilibrate_certified` certifies at 74.1899 cm3 with the
magnesium where GEM-Selektor puts it.

Each missing phase is given what the recipe could make of it,
`min_c b_c / A_cs` over the components it consumes, scaled by a tenth — a
chemical bound, not a guess at the answer. A carbonate in a system holding
1e-9 mol of carbon is offered 1e-9 mol and no more. What matters is only that
the phase starts well away from the boundary, since being *at* the boundary is
what the active-set loop cannot recover from.

### Fixed — a rung is accepted only if it conserves mass

Two things can go wrong at a rung of the continuation, and only one of them
raised. A back end can throw, which was handled; or it can *return* a
composition that violates the element balance, which was accepted and carried
forward as the start of every rung after it. The walk then walks away from the
problem it was posed.

`homotopy_initial_state` now tests each rung before accepting it, and a refused
rung is retaken by halving the distance back to the last `λ` that worked, up to
`max_bisections` times per target — the fixed `steps` ladder becomes a
suggestion, and a rung that cannot be taken in one jump is taken in two.

The test is `|rᵢ| ≤ balance_atol + balance_rtol · scaleᵢ` on the element-balance
residual, row by row, with `scaleᵢ = max(|bᵢ|, Σⱼ |Aᵢⱼ| nⱼ)`, and it **has** to
be mixed rather than relative. Two conservation rows of this problem carry a
legitimately negligible budget: electroneutrality is exactly zero, and a cement
recipe is routinely given a carbon trace of 1e-9 mol. Measured at λ = 0.01 on
the paste, a residual of 1.7e-10 mol on the charge row scores 168 against its
own budget and 1.1e-10 mol on the carbon row scores 10.6 — both physically
nothing. A purely relative criterion rejected 35 of 66 rungs on that walk and
never reached its first three targets; the mixed one accepts all ten rungs and
reaches every target. A row holding 1e-11 mol cannot be balanced better than the
solver's absolute floor, and asking it to be is a category error.

Neither tolerance certifies anything, and both are exposed as keywords rather
than buried as constants. `balance_atol` sits above the accuracy the interior
point itself reaches — about 3e-6 mol on this class of problem — because a rung
is a guess and not an answer. The certificate judges the result, afterwards, and
it is unchanged.

Measured on the CEM I paste at w/c = 0.5, the answer is the same and its balance
is better: certified at 74.1899 cm3 against GEM-Selektor's 74.2136 and
pH 13.0994 against 13.0957, with an element balance of 4.4e-15 where it was
6.3e-14.

The coupled kinetic step is untouched by construction: `implicit_step` passes
`autostart = false`, so it never enters the continuation.

### Fixed — a successful solve printed warnings about its own candidates

A call ending `optimal = true` still printed "equilibrium solve returned
`MaxIters`" and "the dual equilibrium solve did not certify optimality", which
reads as a failed solve and is not one. Those come from **candidates**:
`equilibrate_certified` runs every back end from several compositions precisely
because none of them works on every problem, and keeps whichever answer the
certificate proves, so a candidate falling short is what the search is for.

An internal scope now marks the stretches where starting points are computed or
tried — the back-end solves, the continuation's rungs, the multi-start search
itself, and the coupled kinetic step's warm start — and the two diagnostics stay
quiet inside it. Nothing is hidden: the verdict on the *answer* is still
pronounced once, on its certificate, and `verbose = true` still reports every
rung and every rejected start.

The same distinction fixes something worse than noise. The back-end start solves
are wrapped in a `try`, so under `STRICT_CONVERGENCE[] = true` a `MaxIters`
candidate **raised**, was swallowed, and the search silently lost it — leaving a
caller who asked for strict results with a worse search than one who did not.
Measured on a CEM I paste where the interior point ends on `MaxIters`: with the
flag set the route returned an element balance of 27.6 mol, and with it clear the
very same call returned 1.8e-14. A start is not a result, and the flag now
applies only to results.

### Fixed — the conservation matrix reached the solver with an abstract type

`ChemicalSystem` stores its stoichiometry as `Matrix{Real}` whenever integer and
rational coefficients coexist, which a cement's does — `C3AFS0.84H4.32` and its
kind. That is right for the chemistry and wrong for the solver: an abstract
element type boxes every entry and turns `mul!(res, A, x)` into the generic
fallback with a dynamic dispatch per element, on a product evaluated at every
objective and constraint call. SciMLBase had been saying so on every run of such
a system, warning that "arrays or dicts to store parameters of different types
can hurt performance" — a warning that looked like noise about the library's
internals and was in fact a correct report of a type instability in the hot loop.

`EquilibriumProblem` now narrows the conservation matrix, and the default
`b = A * u0` which inherits the same abstract element type, to a concrete
floating-point array. `DualEquilibriumSolver` had always converted; the
interior-point path had not. A caller who passes a concrete array keeps exactly
what they passed — an exact `Matrix{Rational{Int}}`, a `Matrix{Int}`, or a
`Matrix{<:Dual}` for someone differentiating through it — and the exact
stoichiometry is untouched in `system.SM.A`, where it belongs.

### Fixed — the multi-start search could discard the answer it just computed

`equilibrate_certified` ranks the answers of its multi-start search, and the
comparison was on the **stationarity alone**. A composition can be stationary to
1e-3 while violating mass conservation by moles, and that is not a near-answer:
it is not an answer to the problem posed. Reported from that Windows run, the
first route came back stationary to 2.4e-3 with an element balance off by 6.7 mol
and a phase supersaturated by 45 — and it beat every candidate the continuation
produced, because those were stationary to only 1e-2 while conserving mass. The
search computed a usable answer and threw it away, which is why the certificate
came back **bit-identical** across two successive library fixes.

The ranking is now `_kkt_error` — the worst of stationarity, element balance and
supersaturation, with a negative supersaturation counted as zero since every
absent phase undersaturated is what optimality requires. That is the ranking
`solve_certified` already used internally over its own starts, so the two agree
where they previously disagreed, and there is one definition of the quantity
instead of two.

The certificate of a failed solve now also reports what the automatic initial
approximation did — not reached, declined, produced no usable start, ran without
improving, improved without certifying, or certified. On a solve that fails on
one machine and not another, that is the first thing anyone needs to know, and
it should not require a second run with `verbose = true`.

### Changed — `STRICT_CONVERGENCE[] = true` now raises on an uncertified answer

This is what should have made that Windows run loud instead of plausible. The
flag was read only on the interior-point return code, so a caller who set it —
asking that a non-converged solve never pass as a result — still got an
uncertified answer out of `equilibrate_certified` with nothing but a `@warn`.
An uncertified answer from that route is precisely what the flag exists to
refuse: it can violate mass conservation by moles and still be an ordinary
`ChemicalState`.

**This changes behavior for code that sets the flag**, which is why it is called
out rather than filed under fixes. The default is untouched: with
`STRICT_CONVERGENCE[]` at its default `false`, an uncertified answer is still
returned with a warning, and `optimality_certificate` is still the way to audit
it. Only the opt-in path is affected, and only in the direction the flag asks
for.

One consequence had to be handled inside the package. The coupled kinetic step
computes its warm start with `equilibrate_certified`; that answer is a *starting
point*, not a result, so the flag is cleared around it and restored in a
`finally` — the same distinction the continuation already makes for its own
rungs. Left set, the strict flag would have turned that guess into a raise, the
surrounding `catch` would have swallowed it, and a caller asking for strict
results would have silently got a worse start than a caller who did not.


## v0.15.0 — the alkali end-members of the C-S-H, and a readable aqueous state

The shipped `data/solid_solutions.toml` described a C-S-H that could not hold
alkalis. `CSHQ` was declared with four end-members, and CEMDATA18's `KSiOH` and
`NaSiOH` — the alkali-uptake end-members of that same phase — were reachable
from no file in the package, although both are fully parameterized in both
shipped databases (Cemdata18, Lothenbach et al. 2019; the CSHQ model itself is
Kulik 2011). Uptake of alkalis by the C-S-H is not a refinement: it is what sets
the pore-solution pH of a real paste. `docs/src/examples/cement_wc_ratio.md` said so already, in its
list of limitations — "the alkalis themselves, which in a real paste raise the
pore-solution pH to 13 or above".

Measured against a GEM-Selektor reference on a CEM I at w/c = 0.5 (100 g of
oxides, 50 g of water, the CEMDATA18 phase list, and the extended
Debye-Huckel model that run used): GEMS puts **74 % of the total potassium and
90 % of the total sodium into the C-S-H**. With the six-end-member phase
declared and `HKFActivityModel(Ḃ = 0.097637, Kₙ = 0.0)` on an ion size of zero,
this package now returns pH 13.0994 against GEMS' 13.0957, a total volume of
74.1888 cm3 against 74.2136, `KSiOH` to -0.46 % and `NaSiOH` to +0.06 %, and
activity coefficients within 0.25 % (monovalent) and 1.16 % (divalent) of the
ones GEMS printed. Reaktoro 2.13 on the same problem agrees to 0.08 % on the
volume.

### Added

- **`CSHQ` now has six end-members**: `CSHQ-TobD`, `CSHQ-TobH`, `CSHQ-JenH`,
  `CSHQ-JenD`, `KSiOH`, `NaSiOH`.
- **`C3(AF)S0.84H`**, the CEMDATA18 Fe-siliceous hydrogarnet
  (`C3AFS0.84H4.32` + `C3FS0.84H4.32`), ideal mixing. It is where the iron of a
  ferrite phase ends up, and it took 0.0507 mol in the reference run — second
  only to the C-S-H and the portlandite among the aluminate and ferrite
  hydrates. `build_solid_solutions` therefore returns six phases where it
  returned five.
- The `"C-S-H"` reporting group of `volume_fractions` covers the two alkali
  end-members, which would otherwise have been counted under `"other"`.

### Added — the aqueous state is readable

The activity closures computed the molalities, the ionic strength and the
activity coefficients on their way to the log-activities, and kept all three to
themselves. Anyone comparing a state against GEM-Selektor, PHREEQC or Reaktoro
needs them species by species, and had to reach into
`activity_model(cs, model)` and `ChemistryLab._build_params(state)` to get them.
They are public now, together with the two things that make them meaningful:

- **`molalities(state)`** — `mᵢ = nᵢ / (n_w Mw)` for every solute, mol/kg.
- **`ionic_strength(state)`** — `I = ½ Σ mⱼ zⱼ²`, mol/kg. Model-independent: it
  is a property of the composition, and it is the first thing to compare
  against another code, because an ionic strength that disagrees means the two
  are not describing the same solution whatever their volumes agree on.
- **`activity_coefficients(state, model)`** — γᵢ of every aqueous species,
  evaluated from **the model's own formula** rather than as a ratio `a/m`. The
  ratio agrees for an abundant solute, and the tests check that it does, but a
  species parked at the solver's 1e-16 mol lower bound has its log-activity
  dominated by the closures' `+ ϵ` regularization, and the ratio then returns
  values of order 1e300 for a charge class whose only members are trace. The
  formula depends on the ionic strength and the charge alone, so it is exact at
  any amount. Measured on a CEM I pore solution: the charge classes |z| = 3 and
  4 now come back at -2.7 % and -4.7 % of the coefficients GEM-Selektor
  reports, where the ratio gave 1e300.
- **`log_activities(state, model)`** and **`activities(state, model)`** — the
  vector the Gibbs energy is built from, for every species.
- **`concentration_scale(model)`** — `:molality` or `:molarity`. An activity is
  a number on a scale and nothing in the number says which:
  `DiluteSolutionModel` is on the molarity scale but takes ρ = 1 kg/L, so its
  activities coincide numerically with molalities, while the other two are on
  the molality scale. A custom model should declare its own.
- **`pH(state, model)`** and **`pOH(state, model)`** — `−log₁₀ a(H⁺)` and
  `−log₁₀ a(OH⁻)`, the **activity** convention that GEM-Selektor, PHREEQC and
  Reaktoro report. The existing one-argument `pH(state)` is a different
  quantity: `−log₁₀ c(H⁺)` in mol/L over the computed liquid volume,
  reconstructed through `pKw` when the solution is basic. On a Portland cement
  pore solution at I ≈ 0.2 mol/kg, with γ(H⁺) ≈ 0.61, the two differ by
  **0.21 units** (13.31 against 13.10) — enough to be mistaken for a modeling
  error. Both are now documented side by side.

### Added — the starting point is found, not asked for

`equilibrate_certified` computes an initial approximation when its ordinary
starting points fail to certify, so a realistic cement is solvable without the
caller knowing anything about the answer.

This was the practical obstacle. From the cold state of a CEM I paste of 135
species — all the mass in the reactants, every product at the `ϵ` floor — **no
back end reached the optimum**: the answer came back `optimal = false` with a
worst supersaturation of order 1e1 and a total volume 13 % wrong, and
`equilibrate` exited on `MaxIters` in both `Val(:linear)` and `Val(:log)`. The
problem is convex, so this was never a local minimum; it is the conditioning of
an interior-point method started against the boundary. With 120 of 135 species
at 1e-16 and nine at ~1 mol the barrier gradients span sixteen orders of
magnitude, and a solid-solution end-member has `ln a = ln x → −∞` as its mole
fraction goes to zero, so the objective's gradient is unbounded on exactly the
face where a mixing phase vanishes.

`homotopy_initial_state(state)` walks the solute amount up from a dilute system:
everything but the aqueous solvent is scaled by `λ`, and `λ` goes to 1 with each
step started from the answer to the previous one. At `λ = 1` the composition is
the one given, so the element balance is unchanged; it is a starting point, not
a certified equilibrium. Nothing in it is chemical — "everything but the
solvent" needs no guess at which hydrates will form, and `λ` is not a physical
parameter.

Measured on that paste, `equilibrate_certified(state)` from the cold state now
certifies and returns a total volume of 74.1888 cm3 against GEM-Selektor's
74.2136 (-0.033 %) and pH 13.0994 against 13.0957, identical to what a
chemically informed seed gives. With `autostart = false` the same call fails, at
83.7 cm3.

The route also **restarts from its own answer**. The continuation ends on a
composition that is nearly the equilibrium but not certifiably so, and returning
it as the least-bad answer left the last step to the caller: measured on that
same paste under the per-species Debye-Huckel model, stationarity 9.9e-7 and no
certificate, where one more solve started from that answer gives 1.5e-16 with
the worst absent phase 1.4e-5 below saturation. It is the same observation that
motivates the continuation, applied once more — a start near the answer is what
this problem needs, and the best one available is the answer already in hand.
Bounded, and it stops as soon as a round buys nothing.

Doing that exposed a flaw in how rounds were compared: ranking on the KKT error
whenever the new answer was uncertified let an uncertified point with a smaller
stationarity displace a certified one, trading a proof for a residual. It could
not fire while the incumbent was always uncertified; the restart loop reaches
that rule from a certified state, so the optimality flag is now compared first,
in both directions.

Two design points:

- **It costs nothing in the ordinary case**, because it only runs when nothing
  else certified, and its answer is kept only if it is actually better —
  certified beats uncertified, and among uncertified the smaller KKT error wins.
- **`autostart = false` declines it**, and the coupled kinetic step passes that.
  There the caller already supplies a starting point — the previous instant of
  the integration — and a handful of extra solves inside an implicit ODE step
  would be paid at every step. This is the same principle as a coupled Reaktoro
  run, where the solver is carried across instants.

Differentiability is unaffected: under `ForwardDiff` the certified route strips
to the primal state, solves in `Float64` and attaches the sensitivity through
the implicit function theorem, so the continuation never sees a `Dual` and the
derivative does not depend on how the starting point was found.

The walk is done under `DiluteSolutionModel` whatever the target model is, and
deliberately: walking under the extended Debye-Hückel model with a common ion
size of zero runs away to an ionic strength of 18 mol/kg, its coefficients
falling with `I` raising solubility raising `I`. The ideal model has no such
feedback, and its endpoint is a good start for the non-ideal one.

For the record, the approach this replaced: GEM-Selektor computes its initial
approximation by linear programming — `AutoInitialApproximation` in GEMS3K's
`ipm_simplex.cpp`, an "LPP-based automatic initial approximation of the primal
vector x" using a "modified simplex method with two-side constraints" (Kulik et
al., *Comput. Geosci.* 2013). Reproducing that needs a genuine LP solver. Posing
the same linear program and handing it to the barrier method here does **not**
work: measured, it does not move off the cold state, leaving the linear
objective at -723.8 where -755.9 was feasible. A simplex-based initial
approximation remains the principled option and would need an LP dependency.

### Added — the rest of the CEMDATA18 solid solutions, and a model of any arity

`data/solid_solutions.toml` went from six phases to **eleven**. The five added
complete the set of multi-end-member phases a GEM-Selektor CEMDATA18 run of a
Portland cement is given: `Straetlingite_ss` (straetlingite / straetlingite7),
`AFm_SO4_OH` (C4AH13 / monosulphate12), `AFt_SO4_CO3` (tricarboalu03 /
ettringite03_ss), `Hydrotalcite_AlFe` (Mg3AlC0.5OH / Mg3FeC0.5OH) and `MSH`
(M075SH / M15SH). Every end-member is in both shipped databases.

Ideal mixing is an **assumption** there, and the file says so. CEMDATA18
documents non-ideal parameters for some of these; they are not reproduced
because they could not be sourced with confidence, and a mixing parameter
written from memory is worse than an ideal model honestly labeled. The file also
records what was measured on a CEM I paste: all five are undersaturated
(LogSI −0.032 to −9.981, and a binary ideal solid solution gains at most
log10 2 = 0.301 over its best pure end-member, so mixing cannot bring them in),
and a solid solution whose every end-member sits at the solver's lower bound is
numerically awkward — Reaktoro fails to converge when any single one of the five
is declared on that paste.

**`RegularSolutionModel`** fills a gap in arity, not in chemistry.
`RedlichKisterModel` is the general binary form and is restricted to two
end-members; `IdealSolidSolutionModel` takes any number but no interaction at
all. A C-S-H with six end-members, or the CNASH and ECSH families of CEMDATA18,
had no non-ideal option. The new model is the symmetric multi-component
Margules form, `G^ex = Σ_{i<j} W_ij x_i x_j`, with

    ln γ_k = (1/RT) [ Σ_{j≠k} W_kj x_j − Σ_{i<j} W_ij x_i x_j ]

and `W` in J/mol, symmetric. The tests check the three things that matter: it
reduces **exactly** to `RedlichKisterModel(a0 = W₁₂)` for two end-members, it
satisfies `Σ x_k ln γ_k = G^ex/RT` in a ternary — the Gibbs-Duhem consistency a
hand-written `ln γ` usually gets wrong — and a six-end-member phase accepts it
where Redlich-Kister raises. `build_solid_solutions` reads it from
`model = "regular"` with either `w` (a binary) or `W` (a full matrix).

### Added — two ionic strengths, and a Setschenow coefficient per species

`ionic_strength(state; kind = :effective | :stoichiometric)`. The effective one
(the default, and the previous behavior) sums over the speciated free ions, so a
neutral pair such as `Ca(SO4)@` contributes nothing; it is the one every activity
model here is a function of. The stoichiometric one sums as if every complex were
fully dissociated over the system's primaries, which is the analytical ionic
strength of the recipe and what some salting-out and diffusivity correlations are
fitted against. The gap between them measures how much salt is associated: 0.7 %
on a Portland cement pore solution, far more on a sulfate brine. The
stoichiometric sum reuses the decomposition the mass balance already uses
(`SM.A`), so it needs no separate table of dissociation reactions.

Neither is a difference of *formula* — all three codes compute `½ Σ mⱼ zⱼ²`, and
the spread between GEM-Selektor (0.2097), this package (0.2121) and Reaktoro
(0.2178) on that pore solution comes from their converged compositions, not from
their definitions.

The Setschenow (salting-out) coefficient of a neutral aqueous species is now read
from `sp[:Kₙ]` when present, falling back on the model's global `Kₙ`. CO₂(aq),
the noble gases and the neutral silicates are not equally salted out, and one
coefficient for all of them was the B-dot literature's simplification rather than
a fact. No table of coefficients is shipped: a table is data, and data belongs in
a database or in the caller's hands.

### Verified — the water activity, against Gibbs-Duhem

A benchmark was added for the property that separates this package from its
neighbors. For a 1:1 electrolyte the osmotic coefficient of the extended
Debye-Hückel model has a closed form, and integrating Gibbs-Duhem over the
model's own activity coefficients must reproduce it. `HKFActivityModel` does, and
the test pins NaCl(aq) at 25 °C to a_w = 0.996657, 0.983603 and 0.966898 at 0.1,
0.5 and 1.0 mol/kg — the Gibbs-Duhem-consistent values, which coincide with the
tabulated water activities of NaCl.

Reaktoro's `ActivityModelDebyeHuckel` does not. Measured on the same three
molalities, its reported water activity is 0.20 %, 1.70 % and 3.92 % below what
its **own** activity coefficients imply, so the inconsistency is internal and
needs no external data to demonstrate. On a CEM I pore solution it reports
a_w = 0.8824 where GEM-Selektor gives 0.992588 and this package 0.993723, and the
consequence is quantitative: `ettringite` and `ettringite30` differ by two water
molecules, so their ratio goes as `1/a_w²`, and `(0.9926/0.8824)² = 1.265`
against a measured ratio of ratios of 1.266. That single discrepancy accounts for
the whole difference in the AFt split. Anything whose stoichiometry differs by
water — the C-S-H hydration states, the AFm and AFt series, the hydrogarnets — is
sensitive to it, which is most of a cement.

### Added — a common ion size, settable on the model

`HKFActivityModel(; å)` imposes **one** effective radius on every charged
aqueous species, short-circuiting the per-species tables. That is what
GEM-Selektor, PHREEQC's `-gamma` and most published cement models actually use,
and it was previously reachable only by mutating `sp[:å]` on every species
before building the `ChemicalSystem`. `å_default` looks like the knob for it and
is not: it is the last resort of the lookup chain and is never reached for an
ion either table covers, so setting it changes essentially nothing. Both facts
are now stated in the docstring and asserted in the tests.

`å = 0` collapses the Debye-Hückel denominator to 1, giving the limiting law
plus the B-dot term. With `HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0)`
the package reproduces the activity coefficients of a GEM-Selektor CEMDATA18
run to 0.25 % on the monovalent ions and 1.2 % on the divalent ones — a run
whose model can be recovered from its own output, since CEMDATA18 carries no
ion-size parameter at all.

### Fixed

- `STRICT_CONVERGENCE[] = true` disabled the new initial approximation instead
  of hardening it. The continuation's early rungs are *expected* to fall short —
  they are starting points on the way to `λ = 1`, not results — so under the
  strict flag each one raised, was caught, and every later rung started cold
  again, leaving the walk useless in exactly the mode a careful caller turns on.
  The flag is now saved, cleared for the duration of the walk and restored in a
  `finally`, so a non-converged *result* still raises while an intermediate rung
  does not. Found by running a caller's script with the flag on, and pinned by a
  regression test that checks both the outcome and that the flag comes back.
- `docs/src/tutorials/equilibrium.md` and `README.md` built an AFm solid
  solution from `dict["Ms"]` and `dict["Mc"]`. Neither symbol exists in any
  shipped database, so those examples could never have run. They now use
  `monosulphate12` and `monocarbonate`, which is what the TOML has always used.
- The TOML-format example in `docs/src/tutorials/databases.md` showed the same
  two symbols, and a four-member `CSHQ` that no longer matches the file.
- `docs/src/examples/simplified_clinker_dissolution.md` said "the only built-in
  model is `DiluteSolutionModel`". There are three.

### Documentation

- A new tutorial section,
  [Reading the aqueous properties back](https://micropochemomechanics.github.io/ChemistryLab.jl/stable/tutorials/equilibrium/#sec-aqueous-properties),
  covering the accessors above, the two conventions of pH and why γ is not a
  ratio; and an `Aqueous properties` section in the equilibrium API reference.
- The bibliography gained Helgeson (1969), Helgeson, Kirkham & Flowers (1981),
  Davies (1962), Kulik (2011), Parkhurst & Appelo (2013), Robie & Hemingway
  (1995) and Xu et al. (2011), which the activity-model docstrings had been
  citing in prose without entries. Every DOI was resolved against the Crossref
  REST API; Davies (1962) is a Butterworths monograph that Crossref does not
  index, and its author, title and publisher are confirmed by the
  contemporary review in *Science* (doi:10.1126/science.143.3601.37), while the
  authorship of USGS Bulletin 2131, absent from its Crossref record, is
  confirmed by the USGS publications catalog.
- `data/solid_solutions.toml` now says at its head what it is **not**. `CSHQ`
  and `C3(AF)S0.84H` match the GEM-Selektor phases of those names, but `AFm`,
  `Hydrogarnet` and `Hydrotalcite` are deliberate alternatives to the CEMDATA18
  phase model: GEMS treats `monocarbonate`, `C3AH6`, `C3FH6` and `hydrotalcite`
  as *pure* phases, its AFm solid solution is `C4AH13` + `monosulphate12`, and
  its hydrotalcite solid solution is `Mg3AlC0.5OH` + `Mg3FeC0.5OH` at
  Mg:Al = 3. Reproducing a GEMS result means declaring the phases in the script,
  not taking this file wholesale — and the file no longer lets a reader assume
  otherwise.


## v0.14.2 — holding AMD at a version that still has `SS_Int`

A release with no change to ChemistryLab itself. It exists because an upstream
patch release broke every environment that resolves this package's test target
or its documentation.

AMD 0.5.4, published on 2026-09-04, removed `SS_Int`. `SparseColumnPivotedQR`'s
AMD extension calls it and bounds AMD only at `"0.5.1"`, so any fresh resolve
picks 0.5.4, and `LinearSolve` — hence `OrdinaryDiffEq`, hence this package's
test and docs environments — fails to precompile. Nothing in ChemistryLab is at
fault and nothing in it can avoid the clash.

`AMD = "0.5.1 - 0.5.3"` is therefore declared in the test target and in
`docs/Project.toml`. AMD is **not** a dependency of the package: the bound
constrains nothing downstream, only the environments built here. Remove it once
AMD restores the binding or SparseColumnPivotedQR tightens its own bound.


## v0.14.1 — the ambiguities a test suite cannot see

`Aqua.test_all` now runs as part of the suite. It checks eight properties no
domain test looks at — method ambiguities, unbound type parameters, undefined
exports, dependency hygiene, type piracy, tasks left running at load — and it
found three real defects, all fixed here. No exported name, signature or result
changes; the calls below simply used to fail.

### Fixed

- **Fourteen method ambiguities from one signature.** `convert(::Type{T},
  f::Formula{T}) where {T}` left `T` unconstrained, so dispatch weighed it
  against every `convert(::Type{T}, x)` in the ecosystem — `Missing`,
  `Nothing`, `Ref`, `FunctionWrapper`, `VecElement`, polynomial types and more.
  `Formula{T <: Number}` already implies the bound; saying it in the signature
  costs nothing and removes all fourteen.
- **`convert(::Type{<:ForwardDiff.Dual}, ::Formula)`** was ambiguous with
  ForwardDiff's own catch-all and therefore unusable. It now resolves, spelled
  with the same type parameters ForwardDiff uses so that it is genuinely more
  specific.
- **`promote_rule(Species, Species)`** was the missing diagonal of the
  `Species`/`AbstractSpecies` pair, and errored.
- **`Species(...)` and `CemSpecies(...)`** declared a `where {T}` bound only by
  the varargs, so `T` was unbound when no pair was passed. `T` was never used in
  the body; the bound is gone.
- **`gather_species`** on a dictionary whose keys *and* values are species was
  ambiguous in dispatch and in meaning. It now says which two readings are
  possible instead of raising an unexplained `MethodError`.

### Declared, not fixed

`thermo_factories.jl` extends some thirty `Base` math functions to
`DynamicQuantities.Quantity`, stripping the unit first. That is type piracy —
function and type both belong elsewhere — and it is deliberate. It is now
declared to Aqua as such rather than left unexamined. It remains global: any
code loading ChemistryLab gets these methods whether it asked or not.

## v0.14.0 — an equilibrium that comes with a proof, and a kinetic step that is one problem

### Breaking changes

- `equilibrate(state)` now solves through **every** loaded back end and returns
  the answer `optimality_certificate` proves optimal. It is slower and, on some
  systems, gives a different — correct — answer. `equilibrate(state; certify =
  false)` restores the single-back-end behavior.
- `OptimaSolver = "0.5"` is required. The parameter block the new constraints and
  the kinetic step ride on arrived there; against 0.5's predecessor the extension
  fails with a `MethodError` on `DualNewtonProblem`.
- The `equilibrate` keyword `model` is now named explicitly in the one-argument
  form, and `certify` and `constraint` are new keywords. Code passing positional
  arguments is unaffected.

### An interior-point answer whose charge balance was wrong in the second digit

On 1 mmol of calcite dissolving in 1 kg of water, the interior-point path returned
`‖An − b‖∞ = 3×10⁻⁶`, which looks like round-off. It is not: the rows of a
conservation matrix do not share a scale, and against its own budget that residual
is **3 % on the charge row** and 0.3 % on the carbonate, while the water row at
55.5 mol sits at `1.8×10⁻⁸`. The pH came out 9.5210 against a correct 9.5274.

The measure is fixed in `OptimaSolver` 0.5.0. The answer is fixed here, and the
iteration is not the place: traced over twenty-three iterations, `α` is pinned at
its ceiling of 0.15 and `‖dn‖` decays at `1 − α` with the residual frozen — the
fraction-to-boundary rule lets through 15 % of a correction the next iteration
re-poses. Thirty thousand iterations change nothing, nor does a tolerance of
`10⁻⁶`.

What does is that the problem is convex, so its KKT conditions are sufficient and
`optimality_certificate` **decides** rather than ranks. `equilibrate` therefore
solves by every route and keeps a proved answer:

| case | interior point | dual Newton | certified |
|:--|--:|--:|--:|
| calcite in water, 10–40 °C | 3.0e-2 | 3.0e-12 | proved |
| calcite + 1 mmol CO₂ | 1.0e-3 | 6.3e-13 | proved |
| calcite + 50 mmol CO₂ | 1.3e-16 | 8.5e-14 | proved |
| pure water | 1.3e-16 | 1.3e-16 | proved |
| CEM I, w/c 0.45 and 0.60 | 1.5e-14 | 2.0e-14 | proved |
| CEM I, w/c 0.30 | 7.2e-16 | 3.3e-10, **not** proved | proved |

Ten of ten, against seven through the interior point alone and nine through the
dual Newton alone.

One bug found in the process, and it mattered: `solve_certified` recomputed the
component totals from **each** starting point. A start that violates the balance
therefore posed a different problem, and the dual solve certified the answer to
the shifted one — two "certified" compositions 0.2 % apart on dissolved calcium,
which a convex problem with one minimum cannot have. `b` is now fixed once, from
the state as given.

### Equilibrium under something other than fixed (T, P)

`Adiabatic`, `FixedEnthalpy`, `FixedVolume`, `SealedVolume`, `FixedpH` and
`FixedActivity`, passed as `equilibrate(state; constraint = ...)`.

The prescribed quantity is an unknown of the solver's own square system, not an
outer loop: an adiabatic solve costs one equation, not a solve per trial
temperature. Two vehicles, and they are not interchangeable — a prescribed
**property** adds a parameter and its residual, a prescribed **chemical
potential** adds a *column* to the conservation matrix, an unknown amount of a
titrant the system may draw on, and that amount comes back as part of the answer.

Validated against a published number: adiabatic `H⁺ + OH⁻ → H₂O` gives
**−55.85 kJ/mol** at every amount tested, against the accepted −55.8, with
nothing fitted to it. `ΔT = 6.62 K` for 0.5 mol in a kilogram of water, against
27.9 kJ over 4.18 kJ/K. Enthalpy conserved to `10⁻¹²` relative, and solving at a
*fixed* temperature equal to the one found returns the same composition to
`10⁻¹²`.

A volume constraint on a condensed system is **refused**, with the lever it
measured named in the error. The molar volumes of water and of the minerals in
these databases are exactly pressure-independent, so the relative lever
`(∂V/∂P)·P/V` is about `10⁻⁶` — some 9 600 bar to change the volume by one
percent. That is the physics: the volume of an incompressible condensed system is
fixed by its composition. Declare a gas phase, or use `porosity(state, reference)`
for a sealed specimen at fixed pressure.

### A kinetic step as Leal poses it: one problem, fully implicit

`KineticStepSolver` and `kinetic_step`. The reaction extents are unknowns of the
same Gibbs minimization as the amounts and the element potentials, and the rate is
evaluated at the **end-of-step** composition. No frozen speciation in a right-hand
side, and no lag between the kinetic and the equilibrium species. The reactivity
constraint being linear, it joins the conservation block, so the algebraic cost of
kinetics is the number of **reactions**, not of species.

Measured against an analytic answer: at a constant 1 µmol/s, `Δξ = Δt·r` to
`10⁻¹⁰` relative and the calcite left is `n₀ − kΔt` to the last bit, with the
element balance at `10⁻¹⁴`.

Measured on the feedback path nothing exercised before: with `r = k(1 − Ω)` at
`k = 10⁻⁵ mol/s`, an explicit step of `10⁶ s` would dissolve 10 mol — a thousand
times the calcite present. The implicit step lands at `Ω = 0.999989`, approaching
saturation from below and never crossing it, at any step size.

A **solid solution** takes part: C₃S dissolving into a four-member CSHQ solution
comes out certified at a stationarity of `2×10⁻¹³`, with all four end-members
present and ten times the step giving ten times the C-S-H. Two settings are
decided by the certificate rather than guessed, because neither answer works on
both kinds of problem — `warm_start` for the admission of a mixing phase from a
cold start, `pin_minerals` for whether the kinetic minerals are held in the active
set. Both are documented with the measurements that decide them.

The certificate of a kinetic step is taken on the **augmented** problem. A mineral
held back by a rate law is supersaturated by construction, so tested against the
unconstrained equilibrium it reports a residual of 7 `RT` and a correct answer is
called wrong.

`coupling = :species`, which would impose the assemblage instead of proposing it,
**raises**. It was implemented and does not converge: pinning a pure phase by a
linear row leaves its stationarity row in the system as well, made trivially
satisfiable by that row's own multiplier, and the Jacobian degenerates — the
extents came out 4.3 times short. For a model whose assemblage IS the
stoichiometry, the ODE route already does exactly that, and the tutorial now opens
with a table saying which route answers which question.

### Choosing the step length, and a case the stiff ODE route gets wrong

`kinetic_step_adaptive` chooses `Δt` from a Richardson estimate on the extents —
one step against two half-steps, which for a first-order method is the error of
the coarse one. Reaktoro has no equivalent: the step is the caller's and nothing
reports what taking it cost.

Measured on calcite under `r = k(1 − Ω)`, `k = 10⁻⁴ mol/s`, over `10⁵ s`:

| route | steps | result |
|:--|--:|:--|
| `Rodas5P`, and OrdinaryDiffEq's default polyalgorithm | 6 | **extent −457 mol**, `retcode = Success` |
| `Tsit5`, explicit | 85 626 | correct, 519 s |
| one `kinetic_step` of `10⁵ s` | 1 | wrong, and its certificate says so |
| `kinetic_step_adaptive` | **7** | correct to eight digits |

The cause is not a missing Jacobian term, which is the natural guess. The ODE
route's residual reads the speciation **frozen** at the last accepted step, so
`∂(du)/∂bₑ = 0` is exact for the system being integrated. What goes wrong is that
the frozen speciation makes the right-hand side inconsistent with the state within
a step: an implicit method steps past the point where `Ω` crosses one, the rate
changes sign, and the run enters a branch it never leaves. Three measurements
settle it — bounding `dtmax` to `10³ s`, 103 steps against 6, returns the
identical wrong value to seven digits; removing the re-speciation returns a sane
answer; and routing the re-speciation through the certified route moves −457 to
−383. So the fix is not to freeze, which is what the implicit step does.

The `OrdinaryDiffEq` extension now warns when a trajectory ends on amounts no
chemistry can produce, naming the species and the budget it exceeded. That is
what the ODE route can offer; a rate law that reads the solution belongs on the
implicit one.

### Re-speciation escalates to the certified route when it fails the balance

`respeciate!` went through the interior point alone, and where that is wrong it is
wrong by percent. It now escalates: when the plain solve's absolute element-balance
residual exceeds the retry tolerance, the partition is re-solved through every
back end and the answer the KKT certificate proves optimal is kept. The cost of a
normal run is unchanged, because the escalation fires only on the solves that need
it.

### `ode_solver = :auto`

Hands the choice of integrator to `OrdinaryDiffEq`'s default polyalgorithm, which
detects stiffness at run time. It is offered rather than made the default: on
these problems it lands on the same stiff method as `Rodas5P` and returns the same
answers, so switching would change nothing while removing a reproducible choice.

A step is accepted on the estimate **and** on the certificate. Richardson's
difference is blind to an error the two resolutions share: on calcite over
`10⁵ s` the coarse step and both half-steps each dissolve the entire mineral,
their extents agree to `5×10⁻¹¹`, and the estimator grades that step excellent.
The certificate is not fooled, since such a composition violates
`Δξ − Δt·M·r(n) = 0` by the whole extent. This also closed a hole in
`kkt_certificate`, which checked stationarity, the linear rows and phase
stability but not the **nonlinear** parameter residual, so a step violating its
own rate equation could be proved optimal (`OptimaSolver` 0.5.0).

One further design point, because the obvious choice fails. The tolerance
is relative to the amount each reaction acts on, **not** to the extent. Scaling by
`Δξ` makes the tolerance vanish with the step while the equilibrium solve's noise
does not, so the measured error grows as the step shrinks: the controller halved
to the floor without advancing, and the answer got worse as the tolerance was
tightened — 2.5e-5 at `reltol = 1e-3` against 9.2e-5 at `1e-5`.

### `coupling = :species`: prescribed products inside the strong coupling

The third combination, which neither Reaktoro nor the ODE route offers. One
constraint per kinetic species, `nᵢ = nᵢ(0) + Σⱼ νᵢⱼ Δξⱼ`, so a solid-to-solid
scheme in the form of [Lavergne2018] imposes its products while the aqueous phase
still minimizes `G` under element conservation. Measured on two C₃A pathways, the
extents come out at `Δt·r` exactly and every species follows the stoichiometry to
eight or nine digits, ettringite and monosulfoaluminate included, where
`:reactions` leaves the AFm at zero because the minimization is free to rearrange
the solids.

It was refused in an earlier draft of this release, and finding out why turned up
a genuine defect in `OptimaSolver`: a pinning row for a product that starts absent
has one positive entry and a zero budget, exactly the shape `degenerate_components`
reads as "this component is absent from the system". It pinned that row's
multiplier and declared the species dead — stationarity residual 458, extents 4.3
times short. `conservation_rows` now restricts that criterion to rows that
conserve something.

The pinned species are **eliminated**, not constrained, and that is what makes it
work. Holding them by a linear row — the obvious implementation — stalls: the
species' stationarity row stays in the system, satisfied by a multiplier that must
reach the mineral's own chemical potential, of order 10²–10³ in `RT` units, and
the Newton sits at a fixed point its line search cannot leave. Measured, an
element balance of `6.1×10⁻⁷` mol whatever `maxit`, `tol` or the number of
active-set updates. Their amounts being an explicit affine function of the
extents, they are removed from the system instead, their element content is
subtracted from the budget, and a plain equilibrium is solved over what remains
with a Newton on the `nr` extents around it: the balance goes to `2×10⁻¹⁴` and the
step is certified.

### A single step far beyond the relaxation time can find the other root

For a rate law that vanishes at equilibrium the implicit step has two solutions,
and the second is the composition with the mineral wholly dissolved: it satisfies
the element balance and the reactivity row while violating `Δξ − Δt·M·r(n) = 0`
by the entire extent. The certificate refuses it, and `kinetic_step` reports the
step uncertified rather than passing it off.

Which root the Newton finds is a property of the build, not of the chemistry.
Measured on calcite under `r = k(1 − Ω)` with `k = 10⁻⁵ mol/s`: steps of `10⁴ s`
and `10⁶ s` converge to the right root on Julia 1.12 and to the other one on
1.13.0-rc4 — which is how this was found, a suite green on one and five failures
on the other. `kinetic_step_adaptive` is build-independent, because it refuses an
uncertified step and halves until one certifies; it reaches the equilibrium
values to eight digits in seven steps on either.

The test now asserts the contract rather than the outcome of a Newton on a
particular build: a certified step never crosses saturation, an uncertified one is
flagged and its rate-equation residual is of order one, steps well inside the
relaxation time certify on any build, and the adaptive march reaches the right
answer. Asserting the physics unconditionally on a step of `10⁶ s` was asserting
which root a given LLVM finds.

### `pin_minerals = :auto` ranks the two formulations by margin

It used to return the first one whose certificate passed. That makes the choice a
**threshold** decision, and a threshold decision between two answers that straddle
the tolerance is settled by the last bits: the same commit then behaves differently
on a different machine. Found exactly that way — the suite green locally and five
failures in CI on the identical commit, in the test asserting that a saturating
rate law never crosses equilibrium.

Both formulations are now computed and ranked by a single margin —
stationarity, element balance, worst supersaturation among absent phases, and the
residual of the rate equation. Ranked, the two are not close: on the case that
exposed it the free branch has run to equilibrium and violates its own rate
equation by the whole extent while the pinned one satisfies it to `10⁻¹³`, fifteen
orders of magnitude apart. The adaptive march resolves `:auto` the same way.

### ForwardDiff through the certified route

A composition carrying dual numbers takes the implicit-function route: the
certified solve runs in real arithmetic and the derivative is attached at the
answer. This had to be built rather than inherited — the certified path converts
the component totals to `Float64` to fix them once, which silently dropped every
partial. `∂pH/∂n(CO₂)` now agrees with a central difference to `4.3e-9` relative,
and the primal value is the certified one.

Pushing duals through the search would be wrong in any case: an active set has a
discrete component, so the map `b ↦ n*(b)` is smooth only piecewise and the
derivative belongs at the solution with the active set frozen.

### Correctness of published values

- **The "maleic acid" titration is malonate.** The SLOP98 entries are
  `MALONIC-ACID,AQ` / `H-MALONATE,AQ` / `MALONATE,AQ`, the three-carbon diacid.
  Their `ΔₐG⁰` give pKa **2.851 / 5.696** against a tabulated malonic 2.83 / 5.69,
  and nowhere near maleic's 1.92 / 6.27. The database was right and the name was
  wrong. The script and the page are renamed, and both pKa are now derived rather
  than hard-coded — the old figures were plotted as dashed lines that crossed the
  curve at no half-equivalence, so the defect was visible in the published figure.
- **`parrot_killoh` is deprecated and no longer attributed to Parrott & Killoh.**
  Its nucleation term carries no Avrami logarithm, `K₃` sits where the canonical
  form has `k₂`, `N₁ = 3.3` is the canonical `n₃`, and `k₃ = 1.1` has no
  counterpart: two different models, not two parameterizations. With `PK_PARAMS_*`
  the diffusion branch takes over at α ≈ 0.003 (C₂S), 0.013 (C₃S) and 0.057
  (C₃A), and those three then land on **α(7 d) = 0.2386 whatever their `K₁`**,
  while C₄AF is limited by its own nucleation branch at 0.193. A CEM I at
  w/c = 0.40 is reported near 0.61. All demos and doc pages now use
  `parrot_killoh_avrami` with `PK84_PARAMS_*`: the CEM I paste moves from
  ᾱ(7 d) = 0.234 to **0.628**, ΔT from 2.0 to **14.2 °C**, and the heat released
  from 115 to **308 kJ/kg**. The slag and metakaolin laws move to `waller`.
- **The calorimeter's `Cp` was double-counted.** The denominator is
  `Cp + Σᵢ nᵢ Cp°ᵢ(T)` with the second term recomputed at every step, so adding
  the sample to `Cp` counts it twice and understates ΔT by a factor of about 1.75.
  Three scripts and one doc page corrected; the docstring, which contradicted
  itself, now carries the warning.
- **A documented temperature sweep ran on an empty state.** `ChemicalState(cs)`
  holds no matter, so every element balance was zero and all 21 points returned
  the same trivial answer — flat lines, with prose describing a curve that was not
  there. Repaired, it also shows the note was backwards: `Kₛₚ` is retrograde
  (`10^-8.411` at 10 °C to `10^-8.517` at 30 °C, matching the database to 0.003
  log units, and `−8.48` at 25 °C is the accepted value), and yet dissolved
  calcium **rises**, because the pH falls 0.4 and shifts carbonate to bicarbonate.
  The familiar statement holds for a system buffered at fixed pCO₂, not a sealed
  one.
- **The w/c page was out of navigation, and not by accident.** Its heavy blocks
  were never executed and the scan contradicts its analysis on three points:
  ettringite never forms (the sulfate all goes to `monosulphate12` at
  equilibrium), no clinker survives at any w/c including 0.30, and the porosity has
  neither minimum nor inflection — so no Powers threshold, which is kinetic and
  not thermodynamic. Rewritten with executed blocks, the sealed-curing porosity
  convention (`porosity(eq, fresh)`, 5.3 points and a 7.4 % volume shrinkage apart
  from the one-argument form at w/c = 0.50), its assumptions written out, and put
  back in the navigation.
- Stale documentation removed: the `@test_broken` on `CaOH⁺` at ×2.47 has been
  gone since the convergence test moved to the true KKT error at `μ = 0`, and the
  `pKw = 13.979` that came with it was computed from the pre-fix `OH⁻` — it is
  13.9994 against Reaktoro's 14.0001.

### New tests

`test/published_values.jl` pins the malonate pKa and the retrograde calcite sweep
against the database's own `Kₛₚ`; `test/certified_equilibrium.jl` the certified
route, including that `b` must be fixed once; `test/equilibrium_constraints.jl` the
adiabatic enthalpy of neutralization, the refused volume constraint and the
implicit titrant; `test/kinetics/test_implicit_step.jl` the analytic kinetic step,
the impossibility of crossing saturation, several reactions on one mineral, and a
solid solution inside the step.


## v0.13.0 — a data file is named, not located

### Breaking changes

No name was removed or renamed and no signature changed. Below 1.0 the resolver
treats a minor bump as breaking regardless, so a downstream
`[compat] ChemistryLab = "0.12"` must be widened to `"0.13"`. In this repository
group that is `MeanFieldHomogenization.jl`'s `docs/Project.toml`, which pins
`ChemistryLab = "0.12"` and will otherwise refuse the new version.

One behavior did change, and it can only turn a former failure into a success:
the database readers now fall back to the bundled `data/` directory when a path
does not resolve against the working directory. A call that already worked keeps
resolving to exactly the same file — see below.

### `datapath`, because a script was only runnable from one directory

`build_species("data/slop98-inorganic-thermofun.json")` resolves its argument
against the working directory. That is the repository root when the script is
launched from a shell sitting there, and almost never otherwise: running the
same file from an editor whose REPL started elsewhere failed with

    SystemError: opening file "data/slop98-inorganic-thermofun.json"

Half the scripts in `scripts/` were written that way and half were not; the test
suite had independently converged on `joinpath(pkgdir(ChemistryLab), "data", …)`
in eighteen places, and the documentation pages carried
`"../../../data/…"`, a depth that silently depended on where the page sat in
`docs/src/`.

The exported [`datapath`](@ref) names a bundled file without reference to any
working directory:

```julia
substances = build_species(datapath("cemdata18-thermofun.json"))
ss_phases  = build_solid_solutions(datapath("solid_solutions.toml"), dict)
readdir(datapath())                     # what ships with the package
```

Every script, documentation block and test now uses it, so the code a reader
copies out of a page is code that runs.

### The readers resolve a bundled name, so older code keeps working

`read_thermofun_database`, `build_solid_solutions`, `extract_primary_species` and
`merge_json` (for its two inputs, never its output) resolve their path argument
in this order: as given, then under the bundled `data/` directory, then relative
to the package root. The working directory comes **first**, which is what makes
this safe: a call that resolves today resolves to the very same file tomorrow, a
local database still takes precedence over a bundled one of the same name, and
only a name that *is* one of the shipped files can reach the fallback — so a
mistyped path to a file of your own still fails, now with an error listing what
is available.

The banner these functions print shows the path relative to the package
(`data/cemdata18-thermofun.json`) rather than the absolute one. Documenter
captures that banner into the built page, and an absolute path would have baked
a build machine's directories into the documentation.

### Also

- `scripts/README.md`, which did not exist, records how to run a script and the
  one rule that keeps the collection working: a script consumed by the
  documentation or the test suite (`ionic_hydration.jl`,
  `hydration_calibration.jl`) never calls `Pkg.activate`, because the active
  project is global process state and switching it mid-build breaks every
  `@example` block that follows, on every page.
- The `Usage:` headers of the scripts said `julia --project` while the scripts
  activated `scripts/`; they now agree, and mention that running the file
  directly from an editor works.
- `using Revise` is gone from the three `tutorial_*.jl` scripts — it belongs in
  a `startup.jl`, and it resolved only through the machine's global environment.
  `SymPy`, which `tutorial_symbolic_reactions.jl` imports, is now declared in
  `scripts/Project.toml`, so that script runs in the environment it activates.
- The `merge_json` snippet in `docs/src/man/advanced.md` named a `.dat` file that
  does not exist (`cemdata18.dat`, where the shipped one is
  `CEMDATA18-31-03-2022-phaseVol.dat`) and passed a bundled path as the output
  argument. Both corrected.

### A red test suite, and the version pin that hid it

`test/kinetics/test_calibration.jl` errored with
`UndefVarError: NelderMead not defined in Main`, taking the whole suite down
(496 passed, 1 errored). The cause is upstream: **OptimizationOptimJL 0.4.19
stopped re-exporting Optim's algorithm names**, and
`scripts/hydration_calibration.jl` reached `NelderMead` through that re-export.

What made it look like a local mystery is that the script kept working: a
`scripts/Manifest.toml` is gitignored, and the one on the author's machine still
pinned 0.4.18, where the re-export was there. The test environment resolves fresh
from the registry, gets 0.4.19, and fails — so the suite was red on continuous
integration while every local run of the script was fine.

The script now does `import OptimizationOptimJL: NelderMead`. The binding exists
in that module either way, exported or not, so the explicit import works on both
versions and no longer depends on which one an environment happens to resolve.

## v0.12.0 — calibrated against a measured heat curve, and an isothermal heat that is a heat

### Breaking changes

**The ODE state layout gained one slot for `IsothermalCalorimeter`.** Code that
reads `sol.u` directly, or carries its own offset arithmetic into the state
vector, has to be revisited. `cumulative_heat(sol, ::IsothermalCalorimeter)` and
`heat_flow(sol, ::IsothermalCalorimeter)` also change what they return — from a
reaction extent in moles to a heat in joules. See below.

No name was removed or renamed and no signature changed. Below 1.0 the resolver
treats a minor bump as breaking regardless, so a downstream
`[compat] ChemistryLab = "0.11"` must be widened to `"0.12"`.

### `cumulative_heat` was returning a number of moles

`IsothermalCalorimeter`'s docstring promised that "cumulative heat
`Q(t) = ∫₀ᵗ q̇(τ) dτ` [J]" was "tracked as an extra ODE state". It was not.
`build_u0` appended a trailing slot only for `SemiAdiabaticCalorimeter`, and the
right-hand side wrote `du[end]` under the same condition, so with an isothermal
device `u[end]` was the **last reaction extent ξ_M, in moles**. `cumulative_heat`
returned it as a heat and `heat_flow` finite-differenced it. The function's other
branch, `Q(t) = H(0) − H(t)` from `p.saved_H`, was unreachable: `saved_H` and
`saved_t` are read in exactly two places in the package and written in none.

The isothermal calorimeter now contributes a real state integrating `dQ/dt = q̇`,
so `n_extra_states`, `extend_u0` and `extend_ode!` — until now reachable only from
the test suite — describe what the integrator actually does. The dead branch is
gone. Every index keyed on the state's length was audited; the semi-adiabatic path
already went through `n_extra_states(cal)` and is unaffected, and `p.has_T`
remains false for an isothermal device, so nothing collides on `u[end]`.

A new test compares both routes to the heat on the stoichiometric formulation,
where both are valid, and requires them to agree — the assertion that would have
caught this. The state-layout docstrings were also wrong in a second way, omitting
the reaction-extent block entirely; both now describe all four segments.

One caveat is stated where it belongs rather than left to be discovered. `q̇` here
is `heat_rate`, the heat of the **kinetic** reactions. That is the heat of
hydration when those reactions produce the hydrates. Under partial equilibrium
they only dissolve the anhydrous phases into ions, the hydrates are precipitated
by the Gibbs minimization, and this sum cannot see their heat — worth hundreds of
joules per gram on an ordinary Portland cement. `integrate` now warns on that
combination, as it already did for the semi-adiabatic device, and points at
`heat_release`.

### Calibrating hydration kinetics against measured calorimetry

`scripts/ionic_hydration.jl` has always run a complete CEM I forward with the
*published* Parrott & Killoh parameters, and its `IONIC_CALIBRATION` dictionary
has always been an explicit, deliberately unused hook. Nothing in the repository
closed that loop: there was no experimental dataset, no loss function and no
inverse problem anywhere in `src/`, `scripts/` or `test/`.

`scripts/hydration_calibration.jl` and
[its documentation page](https://micropochemomechanics.github.io/ChemistryLab.jl/stable/examples/hydration_calibration/)
close it, on real measurements, fitting **kinetic parameters only** — the CEMDATA18
thermodynamics are measured quantities, and the Blaine fineness, w/b ratio and
temperature are reported by the experiment.

**The data.** Two CEM I records now ship under `data/experimental/`, subsets of the
CC-BY-4.0 Zenodo deposit of Šmilauer & Reiterman
([10.5281/zenodo.15212785](https://doi.org/10.5281/zenodo.15212785)): a CEM I
52.5 R at w/b 0.50 and Blaine 397 m²/kg over 262 h, fitted, and a CEM I 52.5 R at
w/b 0.45 and Blaine 415 m²/kg over 617 h, held out. They are the only source found
that is at once openly licensed for redistribution, real CEM I measured by
isothermal calorimetry, and documented well enough to *fix* the non-kinetic inputs
— every file states the fineness, the w/b ratio, the temperature and the
normalizing binder mass.

`regenerate.jl` in that directory rebuilds both files **byte for byte** from the
deposit, so nothing about them is unverifiable, and `--check` is an exact
comparison. The subset is the rows nearest to 500 log-spaced times, by nearest
index rather than by interpolation: every number in the vendored files is one the
calorimeter reported. A fourth column carries the depositors' own fitted curve,
read from their tabulated data rather than reimplemented from a functional form
recalled from memory, so a fit can be judged against somebody else's on the same
record. `data/experimental/LICENSE` states the terms — CC-BY-4.0, separate from
the package's LGPL — and the companion article, published CC BY-NC-ND 3.0, is
cited with nothing reproduced from it.

**The published parameters are already close.** With no adjustment at all the
coupled model gives Q(262 h) = 381.9 J/g against 376.0 J/g measured, +1.6 %. The
depositors' four-parameter affinity model, fitted to this very record, gives
388.9 J/g. So parameters fitted in 1984 to other cements, carried through
CEMDATA18 thermodynamics and a Gibbs minimization, land closer at the endpoint
than a purpose-fitted empirical curve. What a calibration can earn here is the
*timing*.

**Three parameters, not six, and the number is measured.** Six candidates were
written down, one per mechanism with a distinct signature on the curve. The
singular values of `∂Q/∂log θ` on the coupled model come out
`[422, 101, 60, 6.3, 1.4, 0.20]` — a factor of nine between the third and the
fourth — so the measurement determines three *combinations*. The correlation
matrix says which: `k₁_C3S` with `n₁_C3S` at −0.985, and `k₃_C3S` with `n₃_C3S`
at −0.977, the latter pair being essentially the whole leading singular direction.
One per pair is all that can be fitted, and the rate constant is kept in each; the
third is `k₁_C3A`. Belite's `k₃_C2S` carries a weight of 0.001 in the three
directions the data see, which is a quantitative way of saying that eleven days of
heat cannot see belite.

Activation energies are not fitted at all, and the reason is not caution: all these
records are at 20 °C, an activation energy is a temperature sensitivity, and a fit
that reported one would be reporting a number the data cannot contain.

**Two failures worth having in the record.** A cheap stand-in for the coupled
model — same rate laws, hydrate assemblage written down by hand — was tried and
does not work, in two distinct ways. First, `C₃A + 3 Gp + 26 H₂O → ettringite` is
stoichiometrically impossible in this mix: 11 % C₃A demands about 1.1 mol of gypsum
per kilogram of binder and 4.6 % gypsum supplies 0.27, so a rate law reading only
C₃A drives the gypsum **negative** and the released heat to 1609 J/g against a
measured 376. Second, with the aluminate sent to hydrogarnet instead, the level is
right to some 9 % but a fit on the *normalized* curve saturates whatever parameter
bounds it is given, whichever subset is fitted — the imposed assemblage produces a
curve shape the Parrott–Killoh family cannot make. Both are now documented sections
of the page and assertions in the test suite, and together they are the
quantitative argument for the model `ionic_hydration.jl` already implements.

**A docstring claim that does not survive measurement.** `parrot_killoh_avrami`
states, after Parrott & Killoh, that "C₃S has no diffusion-controlled stage". The
sensitivity of the released heat to `k₂` for alite was expected to be exactly zero
and is not: the Jander term `k₂(1-ξ)^{2/3}/(1-(1-ξ)^{1/3})` falls as ξ grows, so
past a high degree of hydration it does become the minimum of the three branches,
and `k₂` moves the heat by up to about an eighth of what `k₁` does — all of it
after the first day. The published statement describes the 1984 fit's intent; this
implementation of it has a diffusion-limited tail. Pinned by a test so the claim
cannot drift back.

**Parameter-space automatic differentiation does not work, and the manual said it
did.** `docs/src/man/kinetics.md` claimed "the entire chain is
ForwardDiff-compatible" and showed a derivative through `integrate`. The ODE
right-hand side is AD-clean in the state and in time — deliberately, because
`Rodas5P` needs the time gradient — but not in the parameters: `build_u0` returns
a `Vector{Float64}`; `build_kinetics_params` casts the temperature, the initial
amounts, the stoichiometry and the calorimeter constants to `Float64`; the
surface-area constructors cast to `Float64` despite being declared `{T<:Real}`;
`system_enthalpy` accumulates into a `0.0`; and `respeciate!` is `Float64`-only by
construction. There was no test through `integrate`. The claim is corrected, the
blocking lines are named, and the calibration uses `AutoFiniteDiff()` with
`NelderMead` — still the SciML interface, `Optimization.jl`, and honest about why.

**Reusable, not general.** The helpers live in the script, not in `src/`, on
purpose: the interface should be settled by use before it is exported. The seam for
somebody else's data is `read_calorimetry`, which reads a documented CSV and takes
the experimental conditions from its comment header, and `CALIB_SPEC`, a list of
(phase, rate-law field, prior, box) rows.

### The dormant period is now on by default in `scripts/ionic_hydration.jl`

`run_ionic_hydration` and `ionic_reactions` take `induction` and
`induction_phases`, and the default is **on**: `τ = 5 h`, `m = 2.5`, applied to
the two silicates. A CEM I has a dormant period; Parrot–Killoh does not.

The numbers are round on purpose. The calibration returns 5.6 h and 3.56 on one
record and then shows `τ` correlated with `k₁_C3S` at 0.994, so the data barely
determine either separately. Carrying 20 228 s into an uncalibrated example would
lend a precision the measurement does not support. `induction = nothing` recovers
the previous behavior exactly.

What it changes is the early heat and little else. Measured over 28 days on the
example's own formulation, `Q(6 h)` falls from 74.1 to 40.2 J/g — the point of the
exercise — while the 28-day pore solution (pH 12.754) and porosity (0.3635 against
0.3634) do not move, and the integrator takes 216 accepted steps instead of 200.

Re-measuring both configurations caught an unrelated drift.
`IONIC_DEFAULT_SYSTEM`'s docstring claimed 202 steps and pH 12.58, where the model
*without* any dormant period now gives 200 and 12.754 — that figure had gone stale
for some earlier reason, and the docstring now says so rather than letting the new
default take the blame.

### Downstream: MeanFieldHomogenization.jl

Two things follow for MFH, and the second is not optional.

1. **`docs/Project.toml` pins `ChemistryLab = "0.11"` and must be widened to
   `"0.12"`.** Below 1.0 the resolver treats a minor bump as breaking whatever the
   API did, so the MFH documentation will not resolve against this release until
   that bound moves. MFH's own `Project.toml` does not depend on ChemistryLab at
   all — the coupling lives only in its docs environment.
2. **MFH keeps its own copy of the ionic setup**, in
   `scripts/common/ionic_hydration.jl`, with its own `IONIC_CALIBRATION` and its
   own `parrot_killoh_avrami` call. The dormant period has to be added there, and
   in `scripts/common/stoichiometric_hydration.jl` for the stoichiometric route — on the
   clinker silicates only, since the silica fume there already goes through
   `waller`, whose sigmoid carries its own onset delay.

   It will move a *reported* number, not just an input.
   `scripts/common/paste_micromechanics.jl` returns `E = 0` until the hydrate foam
   percolates — the setting transition as a genuine zero of the self-consistent
   fixed point — and `45_ionic_hydration_micromechanics.jl` reports a setting time
   in hours from the first non-zero modulus. Without a dormant period the model
   accumulates hydrates from `t = 0`, so percolation, and the setting time with it,
   arrives too early.

### `heat_release`'s `states` keyword never worked

`heat_release` accepts `states` so a caller who already holds the certified replay
does not pay for it twice, and its docstring says so: *"it is the expensive part,
and a caller that needs the compositions anyway should not pay for it twice"*. The
line that implemented it was

```julia
states = something(states, speciated_states(sol, kp; times = times))
```

and `something` is an ordinary function, so Julia evaluated **both** arguments and
the replay ran whether or not the states were supplied. Measured on the ordinary
Portland cement of `scripts/ionic_hydration.jl` over 28 days, 60 instants:

| call | before | after |
|--- |--- |--- |
| `heat_release(...; states = sts)` | 121 s | **0.017 s** |
| `heat_release(...)`, replay included | 118 s | 118 s |

`docs/src/examples/ionic_hydration.md` computes its replay once and hands it to
both `ionic_phase_history` and `heat_release` for exactly this reason; it has been
doing the replay twice. A timing assertion in `test/kinetics/test_postprocessing.jl`
now guards it, deliberately — a redundant computation whose result is discarded
cannot be detected any other way.

Note what this does **not** speed up: a calibration loop calling `heat_release`
without `states` already paid one replay and still does. The replay is the price of
correctness, and the same docstring records why — reading the running composition
instead gave 1174 J/g where the certified answer was a few hundred.

Found by benchmarking, not by reading. The hypothesis under test was that the
enthalpy sum dominated; it costs 0.029 s. An optimization written for it was
reverted, because fifty lines and a fallback branch to save 26 ms is not an
improvement, and the measurement that refuted it is what exposed the real defect.

### Proper names out of the identifiers

Author names belong in the bibliography, not in constants. In the two example
scripts and the page that includes one of them:

| was | is |
|--- |--- |
| `LAVERGNE_MIX_C100`, `_VESSEL_CP`, `_LOSS_A`, `_LOSS_B` | `CALORIMETRY_MIX_C100`, `_VESSEL_CP`, `_LOSS_A`, `_LOSS_B` |
| `lavergne_semiadiabatic` | `semiadiabatic_cell` |
| `PK_L2018_C3S` and siblings | `PK_SMOOTHED_C3S` and siblings |

Nothing in `src/` is affected and no exported name changes. `Lavergne2018` and
every citation in prose stay: that is where the attribution belongs.

The names of *models* stay too, and the distinction is deliberate.
`parrot_killoh`, `parrot_killoh_avrami` and `waller` are how the literature
designates those rate laws, in the same way it says Arrhenius, Avrami, Jander,
Powers, Langavant or Blaine — they are technical terms, not credits, and they are
part of the public API. What has been removed is the name attached to a *source of
parameter values*, which is a citation wearing an identifier's clothes.

### Also

  - `ionic_reactions` and `run_ionic_hydration` take a `pk_params` keyword that
    overrides the published Parrott & Killoh sets, so the calibration varies the
    rate laws of the existing model instead of restating it. The default sets are
    now a single `IONIC_PK84` constant rather than a dictionary literal repeated
    per call.
  - `scripts/opc_semiadiabatic_calorimetry.jl` cited Lavergne et al. (2018) with
    the DOI `10.1016/j.cemconres.2017.11.007`, which Crossref resolves to a
    different paper (Machner et al., *CCR* **105**, 1–17). Corrected to
    `10.1016/j.cemconres.2017.10.018`. The same defect was found and fixed in
    `docs/src/refs.bib` in v0.4.0; the script was missed then.
  - Nine bibliography entries added, every DOI resolved against Crossref and the
    Šmilauer entries also against the publisher's own landing page. The dataset and
    the article have different author lists — two creators and three respectively —
    and each is credited for its own artifact.

## v0.11.0 — a proved answer, whichever route found it

### Breaking changes

One new exported name, `solve_certified`. Nothing was removed or renamed and no
existing signature changed, but below 1.0 the resolver treats a minor bump as
breaking regardless, so a downstream `[compat] ChemistryLab = "0.10"` must be
widened to `"0.11"`.

### `solve_certified`: several routes, one proof

`solve_certified(des, starts; b)` solves from each starting composition in turn and
returns the first answer `optimality_certificate` **proves** optimal, or — if none
is proved — the one with the smallest KKT error together with its certificate, so
the caller always sees what it is getting.

Offering several routes is rigorous here rather than opportunistic, and the reason
is the certificate. For a convex problem it does not rank answers, it decides them:
the proof is the same proof whichever start produced the point, and nothing about it
depends on having predicted the winner. What would be a fudge is choosing a route by
taste and reporting its output unproved.

It is also necessary, because no back end dominates. Measured on an LC³ equilibrium
solved COLD at four degrees of reaction, with the dual Newton started from each
interior-point back end in turn:

| degree of reaction | from `OptimaOptimizer` | from `IpoptOptimizer` |
|---|---|---|
| 0.05 | certified, 2.7e-12 | certified, 9.1e-13 |
| 0.25 | **not certified**, 7.2 | certified, 4.2e-12 |
| 0.50 | **not certified** | certified, 4.6e-12 |
| 1.00 | certified, 1.8e-11 | **not certified**, 3.5e-3 |

Each solves a case the other misses. With both offered, every degree of reaction is
certified from a cold start — where before, an intermediate one could only be
reached by continuation from a nearby solution. The starts are supplied by the
caller, so this adds no dependency: whichever back ends are loaded are the ones
available.

### Two assertions corrected because they were wrong

Neither was in the way; both stated something untrue.

The negative control on water autoprotolysis asserted that forcing the Schur
complement back on gives `[H⁺]/[OH⁻] > 2`, pinning a direction of error. Its point
is that the null-space step is what carries the correct ratio, so the assertion is
now that the answer is NOT one — it used to come out at 3.78 and now at 8.7e-5, both
far from unity, and pinning the sign made the test fail for a change that did not
touch what it is about.

The `solve_certified` testset was written asserting `9.5 < pH < 10.3` for a system
that carries 0.01 mol of CO₂. That is the pH of calcite in PURE water, from a
different test; with the CO₂ the answer is 6.43, which is Reaktoro's own value for
this system (`H⁺ = 3.69e-7` in the reference table). The certificate was right and
the assertion was not.

### `speciated_states` walks up to the first instant it is asked for

Every replayed speciation warm-starts from the previous one — and the FIRST one
requested has no previous one, so it started from the cast composition, which
carries no active set at all. The chain is now walked up to it through a few
earlier times whose compositions are discarded and whose only purpose is to hand
the guess an active set.

The consequence was not cosmetic. On the reference OPC replayed at forty instants,
the first of them came back with **56 interior species where the answer has 25** —
every candidate hydrate present, four of them at 1e-5 to 1e-6 mol, the signature of
an interior-point iterate that never reached a vertex — and neither the certifying
Newton nor its continuation recovers from that, because both inherit the start.
That instant is now certified, with stationarity 9.1e-13 and element balance
6.1e-15 mol: the whole replay is proved optimal where before one instant in forty
fell back silently to an uncertified composition. The 28-day heat of that paste
moves from 420.2 to 420.3 J/g, which is the size of the error the fallback was
hiding.

Two other routes to the same instant were implemented and measured, and neither is
in the code: retrying it from each of the thirty-nine certified neighbors, nearest
first, left it unproved with every one of them; and the continuation between
certified instants was rewritten to step forward adaptively rather than bisect —
correct in itself, since bisection lowers the upper end on failure and so abandons
the target for good, but it does not reach this instant either. What was missing was
the start, and only the run-up supplies it.

### The full Portland cement, through its pore solution

New example page and script, `scripts/ionic_hydration.jl`: a complete CEM I —
alite, belite, aluminate, ferrite, gypsum, limestone filler — dissolving into ions,
with the whole hydrate assemblage decided by Gibbs minimization at every accepted
step. The aluminate cascade that a stoichiometric model has to encode by hand comes
out as a result: run with 3.5 % limestone the ettringite survives to 28 days, run
without it the ettringite is depleted to zero, and no line of code distinguishes
the two cases.

The page carries the calorimetry with it, isothermal and semi-adiabatic
(NF EN 196-9), with every number of the cell written out: 420 J/g against 405 J/g
at 28 days, a rise of 19 K against an adiabatic 75 K. The heat is taken from the
enthalpy of the certified compositions, not from the kinetic reactions — those only
dissolve the clinker, so the sum `Σ rᵢ(−Δ_r H⁰ᵢ)` cannot see the precipitation heat
at all, and driving the cell from it gave a rise of 207 K.

## v0.10.0 — the heat of hydration, and solid solutions the solver can hold

### Breaking changes

Five new exported names — `enthalpy`, `heat_capacity`, `missing_enthalpy`,
`system_enthalpy` and `heat_release`. Nothing was removed or renamed and no
existing signature changed, but below 1.0 the resolver treats a minor bump as
breaking regardless, so a downstream `[compat] ChemistryLab = "0.9"` must be
widened to `"0.10"`.

`[compat] OptimaSolver` moves to `"0.4"`, which is required: the solid-solution
work below depends on `SolutionPhase(...; mole_fraction = true)`, introduced
there.

### Calorimetry that a coupled model can actually use

`enthalpy(state)` and `heat_capacity(state)` sum `nᵢ ΔₐH⁰ᵢ(T,P)` and
`nᵢ Cp⁰ᵢ(T,P)` over the composition, and `heat_release(sol, kp; times)` returns
the cumulative heat and the heat rate along a trajectory. This is Eqs. (17)–(21)
of Lavergne et al. (2018): enthalpy is a state function, so its drop between two
states at the same temperature is the heat given off, with reactants, ions and
hydrates each counted once and no reaction stoichiometry to write down.

That last point is what makes it necessary. `heat_rate` sums `rᵢ(−ΔᵣH⁰ᵢ)` over
the KINETIC reactions, which is right when those reactions produce the hydrates —
and wrong under partial equilibrium, where they only dissolve the anhydrous
phases into ions and the hydrates are precipitated by the Gibbs minimization.
Driving a semi-adiabatic cell from it put an ordinary Portland cement at a
temperature rise of 207 K. That combination now warns instead of returning the
number in silence.

`heat_release` reads the **certified** speciations of `speciated_states`, not the
composition the integrator carries. The in-run minimization is warm-started and
uncertified, and a single hydrate is worth hundreds of kilojoules: read that way
the curve came out at 12.7, 145, 1174, 936 and 631 J/g at 1 h, 6 h, 12 h, 1 d and
2 d — heat that rises and then falls. Certified, the same cement gives 185, 288,
347 and 420 J/g at 1, 3, 7 and 28 days, and the curve is monotone.

`missing_enthalpy(state)` lists the species that carry no `ΔₐH⁰` and are
therefore absent from the balance, because a heat curve missing one hydrate is
not visibly wrong.

### Every species now matches Reaktoro, `CaOH+` included

`equilibrium_reference.jl` recorded `CaOH+` as `@test_broken`: out by a factor 2.5
at 4e-9 mol, and left alone because Ipopt landed on the same point, which made it
look like a property of the problem rather than of one solver. It was neither. The
interior-point stage was reporting points it had reached without ever satisfying
the KKT conditions at `μ = 0` — see `OptimaSolver` v0.4.0 — and the trace species
are exactly where that shows. The assertion is now a plain `@test`.

### `s[:Cp⁰]` returned zero on a species that had one

The thermodynamic functions are built on demand. `getproperty` knew that;
`getindex` did not, and returned its not-found value `0` — an `Int64` — for any
of `:Cp⁰`, `:ΔₐH⁰`, `:S⁰`, `:ΔₐG⁰`, `:V⁰` on a species whose functions had not
been forced yet. Callers found out only when they tried to evaluate it. Returning
0 is right for a missing ATOM, which is what that fallback is for; it was never
right for a property the species can produce.

### A solid solution now has a reference worth the name

`DualEquilibriumSolver` took the first end-member of every solid solution as the
phase reference. The reference is the member whose stationarity the outer system
carries rather than inverting, so it has to be the one the phase is mostly made
of — for the aqueous phase that is the solvent, and not by convention. Picking a
minor end-member sends the outer unknown `ln x_ref` towards −∞ and the Jacobian
with it. It is now chosen by magnitude from the caller's own composition, and
solid solutions are declared to `OptimaSolver` as mole-fraction phases, which is
what lets them be solved at all.

## v0.9.0 — equilibria that are proved, not hoped for

### Breaking changes

Two new exported names, `DualEquilibriumSolver` and `optimality_certificate`.
Nothing was removed or renamed and no existing signature changed, but below 1.0
the resolver treats a minor bump as breaking regardless, so a downstream
`[compat] ChemistryLab = "0.8"` must be widened.

### The interior-point solver could not say whether its answer was the answer

`EquilibriumSolver` minimizes `G` by an interior-point method and, on a cement,
it stops on `MaxIters` — at any tolerance. That was known. What was not known is
how much it costs, because nothing in the package could check.

The Gibbs problem is **convex**: an ideal mixing entropy, which is convex, plus
terms that are LINEAR in the amounts of the pure phases, whose potential does not
depend on how much of them there is, over the polyhedron `{A n = b, n ≥ 0}`. The
minimizer is therefore unique and every stationary point is global — so a solver
returning different answers from different starting points is not finding local
minima, it is stopping short of stationarity, and the KKT conditions are
*sufficient*: they can be checked, and checking them is a proof.

**`optimality_certificate`** does that check, on any composition and whatever
produced it. It reports the stationarity of the interior species, the component
balance, and the worst supersaturation among absent phases. A species below the
floor is at its bound, where the condition is the inequality `μ + Aᵀy ≥ 0` and
not the equality — imposing the equality on an amount held at `1e-16` whose
mass-action value is `e⁻³⁰⁰` misstates its log-activity by 263 RT units.

### `DualEquilibriumSolver`

Newton on the KKT system in element-potential space — the Brinkley–Karpov
formulation of the geochemical Gibbs-minimization codes. Writing `u = −Aᵀy`, an
aqueous species obeys `aᵢ = exp(uᵢ − gᵢ)` and a pure phase is present exactly
when `uᵢ = gᵢ` and absent when undersaturated, which is the classical
phase-stability criterion.

**The algorithm is not in this package.** It is
`OptimaSolver.dual_newton_solve`, where it belongs: the active set, the
degeneracy test on the rows of `A`, the two-level Newton and the certificate are
statements about a convex program, not about chemistry. What this package
supplies is the four things that *are* chemistry — the conservation matrix and
the reference potentials `Δ_aG⁰/RT`; the activity model as the callback `h` with
`∇f = g + h(n)`; which species are strictly positive at any solution and which
may vanish; and which species is the solvent, whose activity is a mole fraction
and therefore bounded above, so that its stationarity cannot be inverted. This
release requires OptimaSolver 0.3.

Two levels. The inner one inverts the **solutes'** mass-action laws at fixed
potentials and fixed solvent amount; the outer is a Newton on `1 + m + |P|`
unknowns — the solvent, the `m` element potentials, the amounts of the active
phases — some fourteen numbers for a cement, against forty-seven species in the
interior-point route. Parameterizing the solutes by `ln n` makes their positivity
automatic, so the fraction-to-boundary rule that capped the interior-point step
at *every* iteration has nothing left to act on.

Six things had to be right, and each was found the hard way:

  - **the solvent cannot be inverted through its own mass-action law.** Its
    activity is a mole fraction, so `ln a_w ≤ 0` always, and an arbitrary `y` can
    demand `ln a_w > 0`, for which no finite composition exists. It belongs to the
    outer system, where the balance determines it;
  - **two phases declared saturated over-determine `y`.** Their stationarity rows
    are then jointly infeasible. With portlandite wrongly admitted at 1e-6 mol
    from the starting guess, the residual sat at 29.5 and fell by 1e-3 an
    iteration; dropping it *during* the Newton rather than after, the same solve
    reaches 5e-12 in eleven steps;
  - **an amount must not be clamped to `ϵ` on the way out**, for the reason given
    above;
  - **the active set must be updated during the Newton, not after it**, and one
    phase admitted at a time. Two phases both declared saturated over-determine
    `y` and their stationarity rows are jointly infeasible; admitting a batch
    feeds a cycle in which a phase is admitted, driven negative, dropped, and
    readmitted. On a cement without limestone that produced a solve converged to
    2e-12 — of the wrong subproblem, with an absent phase supersaturated by 10.9.
    Visited active sets are recorded, which makes the loop finite in practice as
    well as in theory;
  - **a component whose total has vanished must be removed from the unknowns.**
    Its element potential is determined by nothing, the Jacobian is singular in
    that direction, and this is the ordinary state early in a hydration run, where
    the iron, sulfur and aluminum totals are at machine zero. The test is not
    `bₖ ≈ 0` but `bₖ ≈ 0` **with the non-zero entries of row `k` sharing a sign** —
    only then does `Σ Aₖᵢnᵢ = 0` with `n ≥ 0` force each term to vanish. The `H⁺`
    row carries `+1` for `H⁺` and `−1` for `OH⁻`, so its zero total is the
    ordinary state of pure water; treating it as degenerate kills the entire
    acid–base system and returns pH 7.000 with the calcite undissolved;
  - **the replay must be seeded and, if need be, continued.** `speciated_states`
    tries a ladder of early instants until one certifies, then walks forward; when
    an instant still resists it bisects the interval in `bₑ` from the last
    certified one, which is a homotopy in the component totals and terminates.
    Without the ladder the first requested instant of a cement could not be
    proved — neither the interior-point answer nor the cast composition puts the
    Newton close enough.

### Measured

On the package's own Reaktoro reference — calcite, CO₂ and water, ideal
activities — the certified answer matches **every species to 1 %**, including
`CaOH+`, which `equilibrium_reference.jl` records as `@test_broken` because the
interior-point answer is 147 % high. On calcite in pure water the certified pH is
**9.90**; the interior-point answer is 6.96, which is not an imprecision but a
wrong answer, and nothing in that solver's output reveals it.

`speciated_states` now certifies each instant it replays, and **names any it
cannot**: falling back silently would leave the caller unable to tell a certified
trajectory from an uncertified one, and those are not the same object. On a full
ordinary Portland cement over 28 days, with and without limestone, **all forty
replayed instants of both mixes are certified**, with element balances between
1e-11 and 1e-13 mol — where the six-hour instant, the AFt peak at which the
assemblage switches, stood at 6.9e-2 mol before this work.

### What is not certifiable yet

**Solid solutions.** `DualEquilibriumSolver` refuses a system that declares one,
and does so explicitly rather than returning a number. A pure phase has unit
activity, hence `gᵢ = uᵢ` and a bound-constrained variable with an active set. An
end-member of a solid solution has `μᵢ = gᵢ + ln aᵢ(n)`, whose activity goes to
`−∞` as its mole fraction goes to zero: it is never exactly absent while the
phase exists, so the active set belongs at the level of the PHASE and not of the
species. That is a different algorithm. Treating an end-member as a pure phase
would not fail loudly — it would return a composition that looks reasonable and
is not the minimum — which is why the refusal is explicit and
[`speciated_states`](@ref) reports such instants as uncertified.

This matters for what comes next: slag and calcined-clay binders form C-A-S-H,
hydrotalcite and AFm solid solutions, and certifying those needs the phase-level
active set.

**The in-run speciation.** `respeciate!` still uses the interior-point solver, so
the compositions *inside* the integration are not certified; only the replay is.
Measured, that makes no difference here, because `bₑ` is integrated from the
rates alone and a Parrot–Killoh or Waller law reads only its own degree of
reaction. A rate law reading log-activities — a saturation ratio, or a
pH-dependent dissolution law — would feed the speciation back into the
trajectory, and would need this closed first.

### Still true

The interior-point solver is unchanged and still does not report convergence on a
cement. It is what gets the certifying Newton into the neighborhood — on 74 of
80 measured cement equilibria the certified answer is reached from its guess —
and that division of labor is deliberate.

## v0.8.2 — a replay that conserves matter, solid solutions that reach the solve

Bug fixes. No API changed and nothing was removed, so a downstream
`[compat] ChemistryLab = "0.8"` accepts this release unchanged.

The `docs/src/examples/coupled_hydration.md` page is rewritten on the real API:
it described four calls that never existed and had no executed block, so nothing
had ever checked it. Every block now runs when the documentation is built.

The `OptimaSolver` bound is relaxed from `"0.2.7"` to `"0.2"`. The patch-level
lower bound guarded against 0.2.6, where `OptimaOptimizer` discarded the
caller's `u0`; a fresh resolution always takes the newest patch, so the bound
only ever mattered when reviving an old manifest.

### The relative residual saturated on near-empty rows

`_row_residual` scaled each row by the larger of its own element total and the
matter flowing through it. For a row whose total is a millionth of the largest —
an element barely present, or the charge row early on — that divisor is
essentially zero, and a rounding error came out as an alarm. A full ordinary
Portland cement balanced to **1.4e-10 mol** at 28 days was reported at `3.2e-2`,
and the near-empty rows of the first steps saturated the measure at `1.0`, which
read as a 100 % violation and was nothing of the sort. Every row is now floored
at a millionth of the system scale, and the same run reports `3.9e-6`.

### A replayed speciation could stop on a broken balance

Two guards were added to [`speciated_states`](@ref), both judged on conservation
of matter rather than on a retcode:

- **A retry from a guess carrying no active set.** The warm start is the previous
  speciation, and it brings the previous *active set* with it. Where the
  assemblage switches — an ordinary Portland cement around six hours, when the
  ettringite peaks and the sulfate starts moving into AFm — that set is the wrong
  one and an interior-point method started inside it does not cross over. When
  the balance exceeds `1e-6 mol` the instant is solved again from the cast
  composition, and whichever result closes better is kept.

- **A feasibility restoration that actually converges.** Alternating projection
  converges linearly, at a rate set by the angle between the box and the affine
  set, and on a cement that angle can be small: at that same six-hour instant,
  where the iron row carries 0.013 mol across thirteen species, 200 sweeps left
  a residual of 8.4e-1 and 20 000 were needed to reach 6.7e-9. A replay runs a
  handful of times and now asks for the budget it needs.

  Inside the ODE right-hand side the budget stays where it was: raising it there
  made the run *worse* — the worst in-run balance went from 1.1 mol to 8.5 at
  2000 sweeps and to 41 at 100 000 — which is the same unpredictability the
  back-end shows elsewhere and is not understood.

Measured on a full OPC replay, with OptimaSolver 0.3: the six-hour instant —
the AFt peak, where the assemblage switches — goes from **6.9e-2 mol to
7.1e-15**, and every other instant lands between 1.8e-15 and 3.6e-9. Over the
201 accepted steps of the coupled run the worst balance is 6.0e-6 mol and the
median 2.7e-9. On the 40-instant grid the chapters use, the porosity and the
elastic modulus are both monotone across every interval and the pore solution
holds at pH 12.52–12.59.

### Solid solutions were dropped from the equilibrium partition

`_equilibrium_subsystem` rebuilds a `ChemicalSystem` for the partition the
equilibrium is actually solved on, and it did not carry `solid_solutions` over.
The loss was silent and total: a coupled run could declare CSHQ, AFm or
Hydrogarnet and the solve would treat their end-members as separate pure phases,
the mixing entropy never entering the Gibbs energy. Measured on alite and belite
with the four CSHQ end-members, the run was bit-identical with and without the
declaration. With the fix the two differ — the element balance goes from 1.5e-12
to 3.6e-15 and the pore solution from pH 12.48 to 12.37.

A solution survives into the partition only if **all** its end-members are in it;
one split across the kinetic and equilibrium sides is not a phase the equilibrium
can mix, and is dropped rather than passed truncated.

This makes the code path work. Whether a given solid solution is *stable* in a
given system is a separate, thermodynamic question: with these Cemdata18
end-members the small silicate system above still precipitates no C-S-H.

### The warning could not tell the trajectory from a probe

The right-hand side is evaluated far more often than the solution advances —
Jacobian finite differences and rejected steps included — and a poor speciation
on a *perturbed* `bₑ` never enters the result. Reporting the worst over all
evaluations alarmed about something that does not affect the answer.

`integrate` now reports the worst on the **accepted steps** first, and the
all-evaluations figure second. It also states when the distinction matters:
`bₑ` is integrated from the rates alone, so a rate law reading only its own
degree of reaction (Parrot–Killoh, Waller) gives a trajectory independent of the
speciation, while a law reading log-activities feeds it back in.

### The warning now gives moles

`integrate` reported only the relative figure, which cannot be judged without
knowing what it was divided by. It now leads with the absolute worst in moles —
`1e-10 mol` is machine precision whatever the system, `1e-2 mol` against a
0.3 mol sulfate budget is not — and says where such cases come from: the first
steps, where the paste has barely reacted. On the OPC above the honest reading
is 1.05 mol at the worst early step and machine precision from three days on.

## v0.8.0 — a public replay of the equilibrium partition

### Breaking changes

- **New exported name: `speciated_states`.** Below 1.0 Julia's resolver treats a
  minor bump as breaking whatever the API did, so a downstream
  `[compat] ChemistryLab = "0.7"` will **not** accept 0.8 and must be widened.
  Nothing was removed or renamed, and no existing behavior changed, so widening
  the bound is the only adjustment required.

### The composition at a given time was not recoverable

`state_at` returns the purely kinetic reconstruction `n(0) + νᵀξ` and says so:
the redistribution performed by the equilibrium solve is not recoverable from
the stoichiometry. For a cement that reconstruction is meaningless — every
dissolved element in solution, not one hydrate. The composition left in the
solver's buffers is no better: it is rewritten at every right-hand-side
evaluation, Jacobian differences and rejected steps included, so it is not the
accepted composition at any time. Reading it is what made a full OPC look as
though it held 0.244 mol of ettringite against 0.267 mol of sulfate.

**`speciated_states(sol, kp; times)`** does the recovery properly: the kinetic
species from the ODE state, and the equilibrium partition re-solved from the
element totals the run carried, walking the instants in order.

Two things it must do, both found the hard way:

- **cap and restore the guess.** The warm start is the equilibrium of the
  *previous* `bₑ`, so once an element has been spent it demands more of it than
  now exists and the interior-point solve begins outside its own feasible set.
- **use a back-end instance the integration has not touched.** Replaying through
  `p.eq_solver.solver` — the object that drove thousands of solves during the
  run — returned pH 14.2 with 0.31 mol of ettringite and no AFm, where a clean
  instance of the same type and settings gives pH 12.58 with the sulfate
  entirely in AFm, on the same trajectory, the same `bₑ` and the same guess.
  **The back-end carries state across solves**; this is a defect upstream, and
  until it is fixed a replay must not inherit it.

The first requested instant is the loose one, having only the cast composition
to start from; every instant after it lands at machine precision. Ask for a
first point within the first hours.

### The back-end defect is fixed upstream

The state the replay had to avoid was `OptimaOptimizer._cache`: with
`warm_start = true` the algorithm object started every solve from its previous
solution, discarding an explicit `u0` — including the guess this package builds
in `respeciate!`, so the caps and projections above were partly defeated during
the run itself. **OptimaSolver 0.2.7** consults the cache only when the caller
supplies no interior point, and `[compat]` now requires it. `speciated_states`
still builds a clean back-end instance, which costs nothing and keeps the replay
correct against an older back-end.

## v0.7.1 — the element balance of a re-speciation, measured and enforced

Bug fixes only. No API was removed or renamed, and no documented behavior
changed, so a downstream `[compat] ChemistryLab = "0.7"` accepts this release
unchanged.

### The reported element-balance residual could not see a violated element

`_respeciate_solve!` scaled the whole residual `|Aₑnₑ − bₑ|` by
`maximum(abs, bₑ)` — in a cement paste, by the water budget. Water is 34 mol
against 0.27 mol of sulfate, so a violation of 0.465 mol of sulfate, i.e. 174 %
of the sulfate present, was reported as `1.4e-2` and read as a converged solve.

The residual is now taken row by row and scaled by each element's own budget,
or by the matter flowing through that row when the budget is zero — as it is for
the charge row, where dividing by `bᵢ` alone turned a rounding error into a
reported residual of 2·10³.

This is a diagnostic change, but not a cosmetic one: every convergence study
run against the old number was measuring the wrong quantity.

### A non-converged speciation was handed on as the next warm start

`eq_warm` was set after any solve that did not throw, non-converged ones
included, so a single bad point seeded every step after it and the error was
locked in for the rest of the run. A speciation is now handed on only when its
per-element residual is within `EQ_RESIDUAL_TOL`; otherwise the step falls back
to the cold stoichiometric reconstruction.

### Re-speciation could start outside its own feasible set

The Gibbs minimization is posed with `Aₑn = bₑ` as a hard equality, and the
warm start is the equilibrium of the *previous* `bₑ`. Once an element has been
spent — the sulfate of an ordinary Portland cement, after the gypsum is gone —
that guess demands more of it than now exists, and the interior-point solve
starts infeasible.

Two guards now run before the solve, both of which leave a feasible guess
untouched:

- `_budget_clip!` caps each species at what the totals can supply,
  `nⱼ ≤ minᵢ bᵢ/Aᵢⱼ`, over the rows it consumes. Rows with a negative total —
  H⁺ in a hydrating cement, met by the hydroxides — bound nothing.
- `_restore_feasibility!` then lands the guess in `{Aₑn = bₑ, n ≥ 0}` by
  alternating projection, ending on the positivity clamp so no small negative
  amount survives for the barrier to reject.

On an isolated ordinary-Portland-cement equilibrium — alite, belite, aluminate
and gypsum dissolved, with ettringite, monosulphate, katoite, portlandite and
C-S-H free to precipitate — this takes the achieved element balance from 174 %
to `2e-15`, machine precision, and returns the textbook assemblage: AFm 0.2323,
ettringite 0.0116, C-S-H 1.8556, portlandite 2.3021, at pH 12.58. Both budgets
close on the last digit — sulfate `0.2323 + 3x0.0116 = 0.2672` and aluminum
`2x(0.2323 + 0.0116) = 0.4879` — against the 0.2672 and 0.4879 available.

### The cold start was supersaturated in every phase at once

The cold-start guess added the stoichiometric reconstruction `Ne^T xi`, which
places every dissolved element in solution with no hydrates at all. For an
aqueous-only system that is harmless; for a cement it is close to the worst
possible start, supersaturated in every phase simultaneously, with an H+ entry
so negative that it is clamped to the floor and loses the acidity the hydroxides
must balance. The guess is now the composition the specimen was cast with, which
`_restore_feasibility!` then carries onto the current `be`.

### What this unblocks

A full ordinary Portland cement now runs end to end: alite, belite, aluminate,
ferrite and gypsum, over 28 days, 202 accepted steps, `retcode = Success`. The
pore solution holds at pH 12.58, and the aluminate sequence comes out of the
thermodynamics with no sequencing rule written anywhere -- ettringite forms
early, peaks at 6 hours and converts to monosulphate once the sulfate is spent,
the AFm settling at 0.26719 mol, exactly the sulfate budget at one SO4 per
formula.

Note that `p.n_full` is a scratch buffer rewritten at every right-hand-side
evaluation, Jacobian differences and rejected steps included; it is not the
accepted composition at `t_end`. Recover a speciation from a solution by taking
`be` from the ODE state and re-equilibrating, as `state_at` already documents.

## v0.7.0 — the real porosity of a setting binder, and a warm-started coupling

### Breaking

- Being a minor bump below 1.0, a downstream `[compat] ChemistryLab = "0.6"`
  will not accept 0.7 and must be widened. No API was removed and no documented
  behavior changed; the additions below are new methods.

### The porosity of a cement was not available

`porosity(state)` returns `(V_liquid + V_gas) / V_total`, and both ends of that
ratio are wrong for a hydrating binder: the denominator is the *current* volume,
which shrinks as the reactions proceed, while a sealed specimen keeps the volume
it was cast with; and the numerator has no gas term, so the empty porosity left
by the chemical shrinkage is structurally invisible. The errors compound — on a
w/c = 0.5 paste at 28 days the method returns 0.327 where the porosity referred
to the specimen is 0.375, the volume having shrunk 7.2 %.

The method is unchanged, and correct for a fixed-volume aqueous system, but now
carries that warning. Two new methods give the right calculation:

- **`porosity(state, reference)`** → `(; liquid, void, total)`, referred to the
  fresh material and counting the Le Chatelier contraction as empty porosity.
- **`saturation(state, reference)`** — a sealed paste desaturates as it hydrates
  though no water ever leaves it.

`porosity(state, ref).void` is exactly the `"void"` entry of `volume_fractions`
called with the same reference, so a transport or micromechanical model fed by
either sees the same thing.

### Fixed — the equilibrium sub-solve was restarted cold at every step

`respeciate!` rebuilt its starting guess from the reaction extents each time,
placing every dissolved element in solution with zero hydrates — a wildly
supersaturated composition. Harmless for an aqueous-only system, nearly the worst
possible start for a cement. It now warm-starts from the speciation left by the
previous accepted step, which is an equilibrium for a nearby `bₑ`.

## v0.6.0 — the coupled path made trustworthy

The v0.5.0 release opened the door to coupled dissolution/precipitation modeling.
Building a cement model on it exposed five defects, four of which were silent.

### Breaking

- **`EquilibriumSolver` gained a `model` field**, so its signature is now
  `EquilibriumSolver{F, S, V, M}`. Code constructing the raw struct positionally
  must be updated; the documented `(cs, model, solver)` constructor is unchanged.
  A new accessor `activity_model(::EquilibriumSolver)` returns it.
- Being a minor bump below 1.0, a downstream `[compat] ChemistryLab = "0.5"`
  will not accept 0.6 and must be widened.

### Fixed — the activity model was silently discarded in a coupled run

A coupled run must rebuild the equilibrium solver for the equilibrium
*sub-system*, because the compiled potential does not carry over. Having no
record of the model it was built from, it rebuilt with `kp.activity_model` —
`DiluteSolutionModel()` by default. Asking for `HKFActivityModel` on a cement
pore solution at I ≈ 0.5 mol/kg therefore gave an infinitely dilute solve, with
no warning. The solver now remembers its model, the rebuild uses it, and a
disagreement with the problem's own model is reported.

### Fixed — the equilibrium sub-solve started on its own bound

`respeciate!` floored its starting guess at `p.ϵ` (1e-30), which
`EquilibriumProblem` raises to exactly the `1e-16` lower bound. An interior-point
method started on its bound stalls: the package's own Reaktoro reference case
reported six non-converged solves out of eight steps for this reason alone. The
guess is now floored strictly inside the box.

**Do not "fix" this class of stall by loosening the optimizer tolerance.**
Measured against Reaktoro 2.13 on that same case, the worst species error grows
from 4.3 % at the default `tol = 1e-10` to 13 % at `1e-9`, 38 % at `1e-8` and
252 % at `1e-7`. A green retcode bought that way costs a factor of sixty in
accuracy.

### Fixed — a bad speciation was invisible

A non-success retcode was a `@warn` at `maxlog = 1` whose result was used anyway,
and it never reached the failure count reported by `integrate`, which only saw
solves that *threw*. Over thousands of steps that is one warning for any number
of bad speciations.

- `NONCONVERGED` counts them, and `integrate` reports the total.
- More usefully, `integrate` also reports the worst `|Aₑn − bₑ|∞` over the run.
  That is the criterion with physical meaning here: a solve can stop short of the
  optimizer's tolerance and still satisfy the element balance to machine
  precision, which is exactly what the remaining stalls do (9e-12 on the
  reference case).

### Fixed — `integrate(kp; reltol = …)` threw a `MethodError`

The documented shortcut forwarded its keywords to a concrete method that accepted
none. Precedence is now explicit: call site beats the solver's own settings,
which beat the defaults.

### Fixed — two solid solutions were dead entries

`data/solid_solutions.toml` declared AFm on `Ms`/`Mc` and Hydrotalcite on
`Ht_OH`/`Ht_CO3`. None of those four symbols exists in either shipped database,
so `build_solid_solutions` skipped both with a warning. AFm being the only
Redlich-Kister entry, the non-ideal mixing path had no live case at all. They are
now `monosulphate12`/`monocarbonate` and, for hydrotalcite, the OH/CO₃ couple at
matching Mg:Al = 2 (`hydrotalcite`/`Mg2AlC0.5OH`). All five phases build.

### Documentation

- `saturation_ratio` stated `ln Ω = Σνᵢlnaᵢ + ln K`, contradicting both the next
  line and the code. The sign is corrected.
- `examples/cement_carbonation.md` claimed `cemdata18-thermofun.json` does not
  contain calcite. It does; the merged database is needed for the phase volumes.
- `man/kinetics.md` used `cs.dict_reactions` as though it were a catalog of the
  database, when it holds only the declared kinetic reactions.

## v0.5.0 — the bridge from a kinetics run to a microstructure

### Breaking

- **The ODE state vector gained the extents of reaction.** Its layout is now
  `[bₑ, nₖ, ξ, (calorimeter)]`, where `ξ` integrates `dξ/dt = r`, one entry per
  kinetic reaction. Code indexing `sol.u` positionally for anything other than
  `bₑ` or `nₖ` must be updated; the calorimeter's slot is addressed from the end
  of the vector and is unaffected, as are `temperature_profile`,
  `cumulative_heat` and `heat_flow`.
- Being a minor bump below 1.0, a downstream `[compat] ChemistryLab = "0.4"`
  will not accept 0.5 and must be widened regardless.

### A rate law may now depend on any species

This is what the state change buys, and it is the point of the release.

Previously only the kinetic species evolved inside the ODE residual: every other
amount was pinned to its initial value for the whole run, because the buffer
holding them was refreshed only by `respeciate!`, which runs only when an
equilibrium solver is present. A rate closure gating on a consumed reactant
therefore never saw it move, and the failure was silent — the kinetic mass
balance stayed exact and the solver reported success while the gated species went
negative. In the cement model that motivated this release, 0.44 mol of ettringite
formed out of 0.27 mol of gypsum.

The residual now reconstructs every species from `n = n(0) + νᵀ ξ` before
evaluating the rates, so reaction sequencing — sulfate available for ettringite,
portlandite available for a pozzolanic reaction, product inhibition — is written
directly. With an equilibrium solver the equilibrium partition is still owned by
the equilibrium solve and refreshed once per accepted step.

As a side effect, [`reaction_extents`](@ref) and [`state_at`](@ref) read the
extents from the solution instead of re-integrating the rates by quadrature: both
are now exact to the solver's tolerance, and the `nsub` keyword is gone.

### From moles to volume fractions

The pieces to turn a hydration run into the input of a mean-field homogenization
scheme were all present — `V⁰` from CEMDATA18, `volume(state, species)`,
`porosity` — but nothing joined them, and the sealed-volume balance was missing
entirely.

- **`volume_fractions(state)`** and **`volume_fractions(state, groups)`**, the
  latter aggregating species into the phase families a homogenization scheme
  consumes (`"C-S-H"`, `"AFt"`, `"anhydrous"`, …). A species listed in two groups
  is an error; species in no group are collected rather than dropped.
- **`chemical_shrinkage(state, state₀)`** and the `"void"` phase. Passing
  `reference` to `volume_fractions` selects the sealed-curing convention: the
  fractions are referred to the initial volume, held fixed, and the deficit left
  by the reactions becomes an explicit gas-filled phase. Without it the fractions
  sum to less than one and the microstructure is wrong.
- **`missing_molar_volumes(state)`** — a species with no `V⁰` is silently absent
  from every volume computation in the package; this makes that visible.

Aqueous solutes have negative standard partial molar volumes, so an individual
fraction can be negative. Those contributions are kept, which is what makes
`volume_fractions` and `volume` agree exactly.

### Post-processing a kinetics solution

`n_full` is a buffer mutated in place — after a run it holds the last accepted
step and nothing else — so there was no way to recover the composition at an
arbitrary time.

- **`reaction_extents(sol, kp)`** — the extent of each reaction, read from the
  ODE state and therefore exact at any instant.
- **`extent_residual`** — the drift between `nₖ` and `νₖᵀ ξ`, which are redundant
  by construction, so the integrator's own consistency is measurable.
- **`state_at(sol, kp, t)`** — the full `ChemicalState` at any instant, every
  species rebuilt from `n = n₀ + νᵀξ` so the result conserves exactly what the
  reactions conserve.
- **`degrees_of_hydration`** and **`mean_degree_of_hydration`**, replacing the
  `phase_alpha` closure copy-pasted into three shipped scripts.

### Parrot & Killoh, canonical formulation

The shipped `parrot_killoh` implements a smoothed variant whose parameters are
not those of the 1984 paper — its diffusion branch uses the same `K₃` and `N₃`
for all four clinker phases. It is unchanged, and now documented as one of two
variants.

- **`parrot_killoh_avrami`** with **`PK84_PARAMS_C3S/C2S/C3A/C4AF`** — the
  canonical Avrami / Jander / power-law form, `α̇ = min(α̇₁, α̇₂, α̇₃)`. With these
  parameters C₂S has no nucleation–growth stage and C₃S no diffusion-controlled
  stage, which is a sharp check on a transcription.
- **`waller`** with **`WALLER_PARAMS_FLY_ASH/SILICA_FUME/SLAG`** — supplementary
  cementitious materials do not follow Parrot & Killoh; `blended_cement_kinetics.jl`
  had to invent PK parameters for slag and metakaolin for want of this.
- **`blaine_factor`**, **`humidity_factor`**, **`powers_alpha_max`** — the three
  rate corrections, previously either absent or retyped inline in every script.

The Avrami branch vanishes at `α = 0`, so `α ≡ 0` solves the ODE and hydration
never starts. Parrot & Killoh's own discrete scheme escapes this by integrating
over the first time step; a continuous solver cannot, so the argument is floored
at `PK_AVRAMI_SEED`.

### Fixed

- **A time-dependent rate law broke every Rosenbrock solver.** The ODE right-hand
  side typed its rate vector from `eltype(u)` alone, but Rosenbrock methods
  (`Rodas5P`, `Rodas4`, `Rosenbrock23`, …) need a *time* gradient, which they take
  by calling the residual with a dual `t` and a plain `u`. Any rate depending on
  `t` then failed with "First call to automatic differentiation for time gradient
  failed". `parrot_killoh` ignores `t`, so nothing exposed it until `waller`. The
  type is now promoted with `typeof(t)`.
- **Non-kinetic amounts were frozen inside the residual** — see the section above,
  which is the substance of this release.

### Bibliography

Every entry of `docs/src/refs.bib` was checked against Crossref (six DOIs) or,
for the four entries that have no DOI, against the publisher or an authoritative
catalog record. Three defects were corrected:

- The DOI of `Lavergne2018` pointed at `10.1016/j.cemconres.2017.11.007`, which
  Crossref resolves to a **different** paper (Machner et al., *CCR* **105**, 1–17).
  Corrected to `10.1016/j.cemconres.2017.10.018`.
- `Lothenbach2015` held the Cemdata18 paper, published in **2019**, while
  `data/solid_solutions.toml` uses the tag `Lothenbach2015` for the genuinely
  different Lothenbach & Nonat (2015) paper and `Lothenbach2019` for Cemdata18.
  The entry is now keyed `Lothenbach2019`, and the real
  Lothenbach & Nonat (2015), *CCR* **78**, 57–70, is added.
- The author of `ParrotKilloh1984` is **Parrott**, with two t's. The citation key
  is unchanged, being referenced throughout the sources and documentation.

## v0.4.0 — Equilibrium actually coupled to kinetics, and dual numbers everywhere

### Breaking

- **`ChemicalState` gained a fourth type parameter**, `R <: Real`, carrying the
  number type of the dimensionless diagnostics: the signature is now
  `ChemicalState{C, S, Q, R}`. Code that spelled the type out in full has to be
  updated; code that used `ChemicalState` without parameters is unaffected.
- `pH`, `pOH`, `porosity` and `saturation` no longer return `Float64`
  unconditionally. They return whatever number type the composition carries, so
  that `∂pH/∂composition` and `∂porosity/∂(w/c)` are obtainable.

### Fixed — the equilibrium was never actually solved during a kinetics run

Four independent defects, each sufficient on its own to disable re-speciation.
They combined to make the whole coupling dead code, and the errors were hidden.

- **`equilibrate(state, ::EquilibriumSolver)` throws.** That is the call
  `build_kinetics_ode` made. `equilibrate` expects a SciML *algorithm* and wraps
  it into an `EquilibriumSolver`; handing it one that is already built wraps it
  twice, and the solve has no applicable method. The entry point for a prebuilt
  solver is `solve(esolver, state)`, which is what is called now.
- **The failure was swallowed by a bare `catch`**, so a run in which
  re-speciation never once succeeded looked exactly like a healthy one.
  Failures are now counted, the first is reported with its exception, and the
  total is reported at the end of `integrate`.
- **`KineticsSolver(; equilibrium_solver = …)` was never read.** `integrate`
  built its parameters from the `KineticsProblem` alone, so the documented way
  of passing a solver did nothing. It is honored now; a solver set on the
  problem still wins, and the conflict is reported.
- **The integrated element amounts `bₑ` were discarded.** The ODE integrates
  `dbₑ/dt = Aₑ νₑᵀ r`, but the re-speciation rebuilt its input from the
  persistent buffer instead, so the constraint of the equilibrium sub-problem
  (Leal et al. 2017, Eq. 54) was not the quantity being integrated. The previous
  speciation is now projected onto `bₑ` through `pinv(Aₑ)` before the solve.

### Fixed — three defects that made the `kinetic_species` API produce nonsense

Found while implementing the partitioned coupling, all on the 3-argument
`KineticsProblem(cs, state, tspan)` path:

- **The controlling mineral was mis-identified.** `_find_mineral_idx` takes the
  first crystalline reactant and, failing that, the first reactant of any kind.
  On a reaction generated from the nullspace — where the kinetic mineral may sit
  on either side — that is the solvent, so the ODE integrated **water** as the
  dissolving phase. A species the system declares kinetic now wins over any
  guess.
- **The stoichiometry had an arbitrary orientation**, so the kinetic mineral
  could carry `ν = +1` and the ODE grew the clinker: α(C₃S) came out at −4.2
  over seven days. Reactions are now normalized to `ν = −1` on their controlling
  mineral, which also delivers the `1/|νₖ|` scaling the `ChemicalSystem`
  docstring already promised.
- **The conservation basis was inconsistent.** `Aₑ` was cut from the canonical
  element matrix `CSM.A` while the solve is posed on `SM.A`, the matrix with
  respect to the primaries; and the equilibrium sub-system re-derived its own
  primaries instead of inheriting the parent's. Either mismatch makes every
  equilibrium solve fail — 89 failures out of 89 steps in one intermediate
  state. `Aₑ` now comes from the sub-system itself, which shares the parent's
  primaries.

With the three fixed, seven days of C₃S/C₂S hydration at `w/c = 0.4` gives
α(C₃S) = 0.277, portlandite and C-S-H in comparable amounts, an alkaline pore
solution at millimolar calcium, and `‖Aₑnₑ − bₑ‖∞ = 8.7e-8`.

### Changed — the coupling is now the partitioned formulation of Leal et al.

`respeciate!` solves `nₑ = φ(bₑ)` over the **equilibrium partition only**, on a
sub-system built for the purpose, with the integrated element amounts `bₑ`
handed to the solver as the constraint rather than derived from a composition.
Running the minimization over the whole system, as the previous code did, would
equilibrate the kinetic minerals instantaneously.

**Validated against Reaktoro.** Calcite dissolving at a constant rate makes the
kinetic half analytic, so every difference of rate-law convention drops out and
only `nₑ(t) = φ(bₑ(t))` is under test. The two codes agree to 0.3 % at `t = 0`
and 4.3 % at 3600 s, the largest deviation always on `CO₂(aq)` at `2.3e-8` mol.
Assertions in `test/coupling_reference.jl`, oracle in
`test/reference/reaktoro_coupling.py`.

### Changed — operator splitting instead of a solve inside the right-hand side

The equilibrium sub-problem is no longer solved inside the ODE right-hand side.
It is solved once per accepted step, by the new `respeciate!`, wired as a
`DiscreteCallback`. The right-hand side is evaluated many times per step and
differentiated to build the Jacobian; solving an optimization problem in there
made the cost unpredictable and the Jacobian inconsistent with the model being
integrated. With the solve outside, residual and Jacobian see the same
speciation.

### Fixed — the conservation matrix was rebuilt by finite differences

`OptimaSolverExt` did not pass `A` and `b` through the problem parameters, so
`OptimaSolver` fell back to reconstructing the constraint Jacobian by finite
differences with a `1e-7` step. The element balance was capped at ~2e-6 mol and
did **not** respond to the requested tolerance. `A` is known exactly; it is
handed over now. Measured on calcite/water: feasibility `2.17e-06 → 4.02e-14`.

This also settles the choice of default back-end. Once the matrix is exact,
OptimaSolver reaches machine-level feasibility *and* is 3 to 26 times faster
than Ipopt on the cement equilibria, so it remains the high-priority default.

### New — differentiating *through* the equilibrium solve

`ForwardDiff` now crosses an `equilibrate`. No solver is asked to iterate on
dual numbers — Ipopt is a C library and never could. The equilibrium is solved
once at the primal values and the sensitivities come from the optimality
conditions of [Leal2017](@cite), the problem Reaktoro solves: one saddle-point
factorization serves every partial derivative, and the answer is exact rather
than a finite difference.

The complementarity block is what partitions the species, and skipping it does
not degrade the answer gently. On calcite + CO₂ in water with a gas phase
declared, the unreduced system puts the whole response into the **absent** gas
species (`n = 5e-11` mol) — a sensitivity that satisfies the element balance to
`4e-16` and means nothing. No back-end returns the stability multipliers `z`, so
the active set is recovered internally: a species negligible on the scale of the
system that nonetheless takes a leading share of the response is pinned and the
system re-solved.

Verified against the package's own finite differences (`9e-5`, the truncation
error) and against **Reaktoro 2.13 reading the same Cemdata18 file**, over the
same eleven species and under the same activity model: the two agree to
`4.4e-4` relative or better on every species, below Reaktoro's own `7.2e-4`
spread across finite-difference step sizes. The equilibrium amounts agree to
the same order. The absent gas species gets exactly zero from the active-set
treatment, against `2e-9` by finite differences.

A cross-code comparison has three knobs — database, species list, activity
model — and each is worth tens of percent. Leaving Reaktoro on HKF against this
package's default `DiluteSolutionModel()` moves `∂Ca²⁺/∂(CO₂)` from `+0.1520`
to `+0.2179`; dropping the aqueous calcium complexes makes `∂Ca²⁺/∂(CO₂)` and
`∂calcite/∂(CO₂)` mirror each other, and restoring them breaks that mirror in
*both* codes alike. The element balance closes to `2e-16` throughout.

### Fixed — trace species now agree with Reaktoro

On calcite + CO₂, every species except `CaOH⁺` agrees with Reaktoro to 5 % or
better, and `pKw` comes out at 13.979 instead of 13.40. Before, `CaOH⁺` was
20.8× high, `OH⁻` 3.19× and `H⁺` 1.24×.

The cause was in the back-end's convergence test, which excluded variables
judged "at their bound" using a threshold scaled by the *largest* amount in the
problem. With a solvent at 55 mol that threshold was `5.5e-7`, so every trace
ion below it was declared to sit on a bound of `1e-16` and its stationarity was
never enforced. `OptimaSolver` 0.2.5 judges it against the variable's own bound.

`CaOH⁺` remains 2.5× high at `4e-9` mol. Ipopt lands on the same value (×2.46),
so it is attributable to neither back-end, and at that amount it is chemically
inconsequential. The assertion is `@test_broken`.

### Changed — the analytic gradient is handed to the back-end

`∂G/∂nᵢ = μᵢ` exactly, the remaining term vanishing by Gibbs–Duhem. The
extension now supplies it instead of leaving the back-end to differentiate
`dot(n, μ(n))`, where that cancellation is only exact in theory: evaluated, the
solvent term alone is `55 × 0.018 ≈ 1` against a `μ` of order 10. This does not
resolve the defect above, but a Gibbs energy is stationary exactly where its
gradient is, and steering on a gradient wrong by ten percent is not defensible
whatever else is true.

### Added — documentation

- **Coupling kinetics and equilibrium** — the partitioned formulation derived,
  including why the state is `(bₑ, nₖ)` rather than `(nₑ, nₖ)`, and what to
  check when a run finishes.
- **The hydrating paste, end to end** — a worked example computing the hydrate
  assemblage instead of imposing it, with measured output.
- **Validation against Reaktoro** — what agrees, what does not, and the three
  knobs a cross-code comparison has to match.

### Fixed — the solver's return code was ignored

Neither extension looked at the optimizer's `retcode`; a non-converged iterate
was written into the state as though it were the equilibrium. It is now checked
and warned about. The warning, rather than an error, is deliberate: the flag is
unreliable in both directions on these problems — pure water comes back flagged
`MaxIters` while giving `[H⁺]/[OH⁻] = 1.000003`. Set
`ChemistryLab.STRICT_CONVERGENCE[] = true` to raise instead.

### Fixed — dual numbers could not enter a `ChemicalState`

The constructor took its element type from the *temperature*
(`Q = typeof(T_q)`), which silently cast the composition to `Float64`: nothing
downstream of a `ChemicalState` was differentiable. The type is now the
promotion of `T`, `P` **and the amounts**. The dimensionless diagnostics follow
(see Breaking above); only their *display* strips the dual part.


## v0.3.1 — Maintenance

- `[compat]` upper bounds raised: `OrdinaryDiffEq` to `"6, 7"`, `TimerOutputs`
  to `"0.5, 1"`.
- CI badge restored; Runic badge.
- Installation instructions updated for registration in Julia's General
  registry; documented the optional optimization backend required to solve
  equilibria (`Optimization`+`OptimizationIpopt`, or `OptimaSolver`).
- Confirmed each GitHub Release keeps archiving automatically to Zenodo's
  existing concept DOI `10.5281/zenodo.17756074` via the native
  GitHub↔Zenodo integration (no workflow or token needed).
- MPCM-Registry deprecated in favor of the General registry.

## v0.3.0 — Chemical kinetics module

### New features

**`src/kinetics/` — mineral dissolution/precipitation kinetics**

- `KineticReaction` — couples a reaction to a rate law; supports
  `transition_state`, `first_order_rate`, and the empirical
  Parrot–Killoh (1984) model for cement clinker hydration
- `KineticsProblem` / `KineticsSolver` — ODE problem formulation
  following the SciML `(u, p, t)` convention; integrated via
  `KineticsOrdinaryDiffEqExt` (weakdep, activated by `using OrdinaryDiffEq`)
- `SemiAdiabaticCalorimeter` / `IsothermalCalorimeter` — coupled
  heat-balance ODE; `temperature_profile`, `heat_flow`,
  `cumulative_heat` post-processing helpers
- `StateView{T,I}` — O(1) named access to species data vectors in the
  ODE hot path (no per-step allocation)
- `RateMechanism`, `RateModelCatalyst`, `BETSurfaceArea`,
  `FixedSurfaceArea` — building blocks for custom rate closures
- `parrot_killoh(params, mineral_name)` factory with built-in
  Schindler & Follliard (2005) Arrhenius correction and default
  parameters for C₃S, C₂S, C₃A, C₄AF
- ForwardDiff-compatible throughout; `KineticFunc` and `transition_state`
  closures accept `Dual` numbers

### Infrastructure
- Relicensed to LGPL-2.1-or-later
- Registered in MPCM-Registry (OptimaSolver resolved via registry)
- GitHub Actions workflows: CI, Documentation, Register, CompatHelper, Format, TagBot
- Multi-version documentation deployment (`docs/deploy_docs.jl`)

## v0.2.3 — Activity models & solid solutions

- Extended Debye–Hückel and Davies activity models
- Redlich–Kister solid solution phases
- `SolidSolutionPhase`, `build_solid_solutions` from TOML
- `with_class` for end-member requalification

## v0.2.0 — Equilibrium solver integration

- `EquilibriumProblem` / `EquilibriumSolver` / `equilibrate` API
- `OptimizationIpoptExt` and `OptimaSolverExt` extensions
- Implicit-differentiation sensitivity via OptimaSolver
- HKF aqueous solute thermodynamic model (`NumericFunc`)
- ThermoFun JSON and PHREEQC `.dat` database parsers
