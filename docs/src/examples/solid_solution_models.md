# [Solid solution models, in numbers](@id sec-app-solid-solutions)

[Solid solutions](@ref sec-theory-solid-solutions) derives the three mixing
models and says what each parameter means. This page evaluates them: what an
end-member's activity actually is, where a regular solution stops being one
phase, and a check that the two non-ideal models agree where they must.

Nothing is solved here either — these are the mixing expressions on their own,
which is the cheapest way to see what a parameter does before putting it in a
`SolidSolutionPhase`.

```@example ss
using ChemistryLab
using Printf

RT25 = 8.31446261815324 * 298.15        # J/mol at 25 °C
nothing # hide
```

## 1. Where a regular solution stops being one phase

A symmetric binary has
``G^{\text{mix}}/RT = x\ln x + (1-x)\ln(1-x) + (W/RT)x(1-x)``, and it unmixes
wherever that is concave. The second derivative at equal fractions is
``4 - 2W/RT``, so the critical point sits at ``W = 2RT`` — about 5 kJ/mol at
25 °C:

```@example ss
d2G(x, w) = 1 / x + 1 / (1 - x) - 2w          # w = W/RT
verdict(v) = abs(v) < 1.0e-12 ? "critical point" : v > 0 ? "one phase" : "unmixes"
println("  W/RT    W (kJ/mol)    ∂²(G/RT)/∂x² at x = 0.5    verdict")
for w in (0.0, 1.0, 1.9, 2.0, 2.1, 3.0)
    @printf("%6.2f   %10.2f   %22.3f    %s\n",
            w, w * RT25 / 1000, d2G(0.5, w), verdict(d2G(0.5, w)))
end
```

This has to be checked before a parameter goes into a model, because nothing in
the equilibrium solver detects it: because a solid solution is entered as **one**
phase and the activity expression keeps returning numbers on the wrong side of
the threshold. Those numbers describe a metastable single phase.

## 2. Redlich-Kister reduces to a regular solution, as it must

``a_0`` is the symmetric term of the Redlich-Kister expansion, so it is exactly
a regular solution's ``W_{12}``, and setting the asymmetric terms to zero must
reproduce the other model. Two independently written methods, one identity:

```@example ss
W = 12_000.0        # J/mol
T = 298.15
x = [0.3, 0.7]

reg = RegularSolutionModel([0.0 W; W 0.0])
rk = RedlichKisterModel(a0 = W, a1 = 0.0, a2 = 0.0)

for k in 1:2
    lr = ChemistryLab._excess_ln_gamma(reg, k, x, T)
    lk = ChemistryLab._excess_ln_gamma(rk, k, x, T)
    @printf("end-member %d:  regular ln γ = %+.10f   Redlich-Kister ln γ = %+.10f   Δ = %.2e\n",
            k, lr, lk, abs(lr - lk))
end
```

This ``W`` corresponds to ``W/RT = 4.8``, past the critical value of 2 from §1.
The identity is algebraic and holds regardless of stability, but the composition
it is evaluated at is not one a homogeneous phase would occupy.

## 3. What the models do to an activity

``W = \pm 4`` kJ/mol below, which is ``W/RT = 1.6`` — inside the stable range on
purpose, so the numbers describe a phase that stays homogeneous:

```@example ss
models = ["ideal" => IdealSolidSolutionModel(),
          "regular W>0" => RegularSolutionModel([0.0 4000.0; 4000.0 0.0]),
          "regular W<0" => RegularSolutionModel([0.0 -4000.0; -4000.0 0.0])]

@printf("W/RT = %+.2f, so §1 says one phase everywhere\n\n", 4000.0 / RT25)
println("            a₁ = x₁ γ₁")
println("   x₁      ideal    W>0      W<0")
for x1 in (0.01, 0.1, 0.3, 0.5, 0.7, 0.9)
    xs = [x1, 1 - x1]
    a = [x1 * exp(ChemistryLab._excess_ln_gamma(mod, 1, xs, 298.15)) for (_, mod) in models]
    @printf("%6.2f   %7.4f  %7.4f  %7.4f\n", x1, a...)
end
```

Two consequences. A positive ``W`` pushes the activity of a dilute end-member
**above** its mole fraction — the host is rejecting it — while a negative ``W``
pulls it below, the host stabilizing it. And at ``x_1 \to 1`` all three converge,
because ``\gamma_1 \to 1`` as the phase becomes pure: the standard state of an
end-member is the pure end-member, and the models agree there by construction.

The ideal column is the relevant one for cement. An end-member at
``x = 10^{-3}`` has ``a = 10^{-3}``, so its saturation index is shifted three
decades below the pure phase — a trace component is stabilized simply by being
diluted in a host, which is how a solid solution takes up an ion that would
never precipitate on its own.

See also: [Solid solutions](@ref sec-theory-solid-solutions) for the derivations,
and [What the choice of activity model costs](@ref sec-app-activity-models) for
the same exercise on the aqueous side.
