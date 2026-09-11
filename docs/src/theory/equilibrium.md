# [Proving that an answer is the answer](@id sec-theory-certificate)

An interior-point method minimizes ``G`` by walking the interior of the feasible
set, and on a cement equilibrium it stops on `MaxIters` — at any tolerance.
Whether the point it returns is the minimum is then an open question. This page
is why the question has an answer at all, what a certificate checks, and how a
solver can aim at the optimality conditions directly instead of at the objective.

The formulation being certified is the one set out in
[Thermochemistry](@ref sec-theory-thermo) §4: minimize ``\sum_i n_i g_i(n)``
subject to ``\mathbf{A}n = b`` and ``n \ge 0``.

## Why the question has an answer

Write the Gibbs energy in `RT` units as `G(n) = Σᵢ nᵢ μᵢ(n)`. Its ideal part

```math
\varphi(n) = \sum_i n_i \ln\frac{n_i}{N}, \qquad N = \sum_j n_j
```

has Hessian ``\operatorname{diag}(1/n_i) - \tfrac1N \mathbf{1}\mathbf{1}^\top``,
and for any ``v``

```math
v^\top \nabla^2\varphi\, v = \sum_i \frac{v_i^2}{n_i} - \frac{1}{N}\Bigl(\sum_i v_i\Bigr)^2 \;\ge\; 0
```

by Cauchy–Schwarz applied to ``v_i = (v_i/\sqrt{n_i})\sqrt{n_i}``. A pure phase
has unit activity, so it contributes a term **linear** in its amount. Hence `G`
is convex, the feasible set ``\{An=b,\ n\ge 0\}`` is a polyhedron, and — the
constraints being affine, so that the linearity constraint qualification holds
everywhere — the KKT conditions are **necessary and sufficient**.

Two consequences follow. The minimizer is unique, so a solver returning different
answers from different starting points is not finding local minima but stopping
short of stationarity. And optimality can be *checked*: a composition satisfying
the KKT conditions is proved globally optimal.

## The certificate

[`optimality_certificate`](@ref) checks the three conditions, on any composition
and whatever produced it. Writing ``u = -A^\top y`` for the element potentials:

| condition | on which species | meaning |
|:--|:--|:--|
| ``\mu_i + (A^\top y)_i = 0`` | interior (`n > floor`) | stationarity |
| ``An = b`` | — | conservation of matter |
| ``u_i \le g_i`` | a **pure** phase at its bound | that phase undersaturated |
| ``\ln \sum_i \exp(u_i - g_i - \ln\gamma_i) \le 0`` | a **mixing** phase held entirely absent | that solution cannot form |

The last row is Michelsen's tangent-plane measure, and it is a separate test
because a mixing phase needs one: its members are never exactly zero while it
exists, so they are neither interior nor at a bound, and a solid solution left out
of the assemblage used to pass the certificate **unexamined**. The trial
composition is refined against the phase's own activity model, so the test is not
the ideal approximation.

Two subtleties decide whether the check is meaningful.

A species **at its bound** obeys the inequality, not the equality. Imposing the
equality on an amount held at `1e-16` whose mass-action value is `e⁻³⁰⁰`
misstates its log-activity by 263 `RT` units, and the check then reports a
residual of 74 for a composition solved to `5e-12`.

A species carrying a **vanished component** is absent by the *constraint*, not by
thermodynamics, and its saturation index is meaningless — the element potential
of a component nobody supplies is determined by nothing. The test for that is not
`bₖ ≈ 0` but `bₖ ≈ 0` **with the non-zero entries of row `k` sharing a sign**:
only then does ``\sum_i A_{ki} n_i = 0`` with ``n \ge 0`` force each term to
vanish. The `H⁺` row carries `+1` for `H⁺` and `−1` for `OH⁻`, so its zero total
is the ordinary state of pure water; treating it as degenerate kills the entire
acid–base system and returns pH 7.000 with the calcite undissolved.

## The certifying solver

[`DualEquilibriumSolver`](@ref) solves the KKT system directly, in element
potentials. From ``\mu_i + (A^\top y)_i = 0`` an aqueous species obeys the
mass-action law ``a_i = \exp(u_i - g_i)``, and a pure phase is present exactly
when ``u_i = g_i``, absent when undersaturated — the classical phase-stability
criterion.

Two levels. The inner one inverts the **solutes'** mass-action laws at fixed
potentials and fixed solvent amount; the outer is a Newton on ``1 + m + |P|``
unknowns — the solvent, the `m` element potentials, and the amounts of the
active phases. Parameterizing the solutes by ``\ln n`` makes their positivity
automatic, which is what removes the fraction-to-boundary limit that caps the
interior-point step at every iteration.

The solvent is deliberately **not** inverted through its own mass-action law: its
activity is a mole fraction, so ``\ln a_w \le 0`` always, and an arbitrary `y` can
demand more, for which no finite composition exists. It belongs to the outer
system, where the balance determines it.

!!! note "What it buys, measured"
    On calcite in pure water the certified pH is **9.90** against an
    interior-point 6.96 — not an imprecision but a wrong answer, and one nothing
    in that solver's output reveals. On the Reaktoro reference (calcite, CO₂ and
    water) both routes now agree with Reaktoro on every species: above `10⁻⁵` mol
    to `10⁻³` relative, the trace ions to 5 %, the worst being `CaOH⁺` at ×1.032
    on 1.6 nmol. That reference used to carry a `@test_broken` for `CaOH⁺` at
    ×2.47; what closed it was the convergence test moving to the true KKT error at
    `μ = 0` (`OptimaSolver` 0.4.1), and `test/equilibrium_reference.jl` is now 26
    plain assertions.

    [`speciated_states`](@ref) certifies every instant it replays and names any it
    cannot. On a full ordinary Portland cement over 28 days, all forty replayed
    instants are certified, with element balances between `1e-11` and `1e-13` mol.

## See also

  - [Thermochemistry](@ref sec-theory-thermo) — the minimization and its
    multipliers
  - [Chemical Equilibrium](@ref sec-equilibrium) — driving the certifying route
  - [`equilibrate_certified`](@ref), [`optimality_certificate`](@ref),
    [`DualEquilibriumSolver`](@ref), [`saturation_indices`](@ref)
