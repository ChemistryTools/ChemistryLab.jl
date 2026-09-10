# [Solid solutions](@id sec-theory-solid-solutions)

A pure solid has ``a = 1``: it is its own standard state, its chemical potential
does not depend on how much of it there is, and it either precipitates or it does
not. A **solid solution** is a single crystalline phase whose composition varies
continuously between end-members, and that changes the thermodynamics
qualitatively: the activity of an end-member depends on the composition of the
phase, so the phase can absorb a component gradually instead of appearing all at
once.

This matters in cement more than anywhere else. C-S-H is not a compound with a
formula but a phase of variable Ca/Si; the AFm phases exchange sulfate, carbonate
and hydroxide on the same interlayer site; alkalis partition into hydrates rather
than staying in solution. A model that only knows pure phases predicts sharp
appearances and disappearances that are not observed.

```@example ss
using ChemistryLab, Printf
nothing # hide
```

## 1. Ideal mixing: the activity *is* the mole fraction

Put ``n_k`` moles of each end-member into one phase, with mole fractions
``x_k = n_k/\sum_j n_j``. Even with no interaction at all, mixing lowers the
Gibbs energy, because there are more ways to arrange a mixture than a pure
solid. Per mole of phase,

```math
G^{\text{mix}} = RT\sum_k x_k \ln x_k \;\; (<0) ,
```

which is entropy alone — the ``-T\Delta S`` of mixing. Differentiating for one
end-member gives its chemical potential and hence its activity:

```math
\mu_k = \mu_k^\circ + RT\ln x_k
\qquad\Longrightarrow\qquad
\boxed{\;\ln a_k = \ln x_k\;}
```

That is [`IdealSolidSolutionModel`](@ref), and it is the default. Note what it
already buys: an end-member at ``x_k = 10^{-3}`` has ``\ln a_k = -6.9``, so its
saturation index is shifted by ``-3`` decades relative to the pure phase. A trace
component is *stabilized* by being diluted in a host, which is why a solid
solution can take up an ion that would never precipitate as its own phase.

## 2. Non-ideal mixing: the excess Gibbs energy

Real end-members interact. Everything beyond the ideal term is collected into
the **excess** Gibbs energy ``G^{\text{ex}}``, and the activity coefficient is
its partial molar derivative:

```math
G = \underbrace{\sum_k x_k\mu_k^\circ}_{\text{mechanical mixture}}
  + \underbrace{RT\sum_k x_k\ln x_k}_{\text{ideal mixing}}
  + \underbrace{G^{\text{ex}}(x,T)}_{\text{interaction}} ,
\qquad
\ln\gamma_k = \frac{\partial}{\partial n_k}\!\left(\frac{n\,G^{\text{ex}}}{RT}\right) ,
```

and then ``\ln a_k = \ln x_k + \ln\gamma_k``. This is exactly what the code
does: `_solid_solution_lna!` writes `log(x[k] + ϵ) + _excess_ln_gamma(model, k, x, T)`
into the log-activity vector, so the two terms are visibly separate and a model
supplies only the second.

## 3. `RegularSolutionModel` — one interaction energy per pair

The simplest non-ideal form keeps only pairwise contacts, with one energy per
pair of end-members:

```math
\frac{G^{\text{ex}}}{RT} = \sum_{i<j}\frac{W_{ij}}{RT}\,x_i x_j ,
\qquad
\ln\gamma_k = \sum_{j\neq k}\frac{W_{kj}}{RT}x_j \;-\; \sum_{i<j}\frac{W_{ij}}{RT}x_i x_j ,
```

which is the symmetric multicomponent Margules expression, and it is
`_excess_ln_gamma(::RegularSolutionModel, …)` line for line. For a **binary** it
collapses to the form usually quoted:

```math
\ln\gamma_1 = \frac{W_{12}}{RT}x_2^2 ,
\qquad
\ln\gamma_2 = \frac{W_{12}}{RT}x_1^2 .
```

``W_{ij}`` (J/mol) is the energy of an ``i``–``j`` contact *relative to the
average* of ``i``–``i`` and ``j``–``j``, so its sign is physical:

  - ``W_{ij} > 0`` — unlike neighbors are unfavorable. ``\gamma > 1``, the phase
    resists mixing, and above a threshold it unmixes;
  - ``W_{ij} < 0`` — unlike neighbors are favorable. ``\gamma < 1``, mixing is
    stabilized beyond ideal, and the phase may order.

`W` must be square and symmetric, which a **validating inner constructor**
enforces; the reason it is inner rather than outer is recorded in the source, and
it is a mistake worth not repeating.

### The threshold is ``W = 2RT``, and it is checkable

For a symmetric binary, the molar Gibbs energy of mixing is
``G^{\text{mix}}/RT = x\ln x + (1-x)\ln(1-x) + (W/RT)x(1-x)``. A phase is
unstable to unmixing where that function is concave, and

```math
\frac{\partial^2}{\partial x^2}\frac{G^{\text{mix}}}{RT}
 = \frac{1}{x} + \frac{1}{1-x} - \frac{2W}{RT} ,
```

which at ``x = 1/2`` is ``4 - 2W/RT``: negative as soon as ``W > 2RT``. So
``W/RT = 2`` is the critical point of a symmetric regular solution — below it one
homogeneous phase, above it a miscibility gap that widens as ``W`` grows:

```@example ss
d2G(x, w) = 1 / x + 1 / (1 - x) - 2w          # w = W/RT
verdict(v) = abs(v) < 1.0e-12 ? "critical point" : v > 0 ? "one phase" : "unmixes"
println("  W/RT     ∂²(G/RT)/∂x² at x = 0.5     verdict")
for w in (0.0, 1.0, 1.9, 2.0, 2.1, 3.0)
    @printf("%6.2f   %22.3f     %s\n", w, d2G(0.5, w), verdict(d2G(0.5, w)))
end
```

## 4. `RedlichKisterModel` — a binary that need not be symmetric

A regular solution is symmetric by construction: swapping the end-members leaves
``G^{\text{ex}}`` unchanged. Real binaries often are not, so the excess is
expanded in a polynomial of the composition difference. As implemented, for two
end-members,

```math
\ln\gamma_1 = \frac{x_2^2}{RT}\Big[a_0 + a_1(3x_1 - x_2) + a_2(x_1 - x_2)(5x_1 - x_2)\Big] ,
```
```math
\ln\gamma_2 = \frac{x_1^2}{RT}\Big[a_0 - a_1(3x_2 - x_1) + a_2(x_2 - x_1)(5x_2 - x_1)\Big] ,
```

with ``a_0, a_1, a_2`` in J/mol. Read the roles off the formulas: ``a_0`` is the
symmetric term and is exactly a regular solution's ``W_{12}``; ``a_1`` and
``a_2`` are the asymmetric corrections, and setting them to zero must reproduce
the regular model. That is an identity between two independently written methods,
so it is worth testing rather than trusting:

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

`RedlichKisterModel` **requires exactly two end-members**, and
[`SolidSolutionPhase`](@ref) refuses the combination at construction rather than
letting a ternary reach an expression written for a binary. For three or more
end-members the choices are the ideal model or `RegularSolutionModel` with a
full ``W`` matrix.

## 5. What the three models do to an activity

Same binary, same ``W``, the whole composition range:

``W = \pm 4`` kJ/mol is used below, which is ``W/RT = 1.6`` at 25 °C — inside
the stable range of §3 on purpose, so that the numbers describe a phase that
actually stays homogeneous:

```@example ss
RT25 = 8.31446261815324 * 298.15
models = ["ideal" => IdealSolidSolutionModel(),
          "regular W>0" => RegularSolutionModel([0.0 4000.0; 4000.0 0.0]),
          "regular W<0" => RegularSolutionModel([0.0 -4000.0; -4000.0 0.0])]

@printf("W/RT = %+.2f, so §3 says one phase everywhere\n\n", 4000.0 / RT25)
println("            a₁ = x₁ γ₁")
println("   x₁      ideal    W>0      W<0")
for x1 in (0.01, 0.1, 0.3, 0.5, 0.7, 0.9)
    xs = [x1, 1 - x1]
    a = [x1 * exp(ChemistryLab._excess_ln_gamma(mod, 1, xs, 298.15)) for (_, mod) in models]
    @printf("%6.2f   %7.4f  %7.4f  %7.4f\n", x1, a...)
end
```

Two readings. A positive ``W`` pushes the activity of a dilute end-member
*above* its mole fraction — the host is rejecting it — while a negative ``W``
pulls it below, which is the host stabilizing it. And at ``x_1 \to 1`` all three
converge, because ``\gamma_1 \to 1`` as the phase becomes pure: the standard
state of an end-member is the pure end-member, and the models agree there by
construction.

!!! warning "A large positive `W` can put the phase outside its own stability"
    With ``W = 12`` kJ/mol — the value used for the identity check in §4 — the
    ratio is ``W/RT = 4.8``, well past the critical 2, so a *homogeneous* binary
    at ``x_1 = 0.5`` is not what that model describes: it would unmix. The
    activity expression keeps returning numbers there, and they are the numbers
    of a metastable single phase. Nothing in the equilibrium solver detects it,
    because a solid solution is entered as one phase; checking ``W/RT`` against
    §3 is the caller's job.

## 6. How a solid solution is declared

A [`SolidSolutionPhase`](@ref) names its end-members and carries a model:

```julia
phase = SolidSolutionPhase("CSHQ", [sp["CSHQ-JenD"], sp["CSHQ-JenH"],
                                    sp["CSHQ-TobD"], sp["CSHQ-TobH"]];
                           model = IdealSolidSolutionModel())
```

Every end-member must be `AS_CRYSTAL`, and the phase is checked against its
model's arity at construction. Sets of end-members with their interaction
parameters are read from a TOML through
[`build_solid_solutions`](@ref), which is also the only interaction-parameter
file format in the package — `data/solid_solutions.toml` is the shipped example.

Two consequences to keep in mind when reading results:

  - a `ChemicalSystem` groups end-members into `ss_groups`, and **every** aqueous
    activity model must call the solid-solution branch or the end-members come
    back with ``\ln a = 0``, i.e. treated as pure phases;
  - an end-member's [`saturation_indices`](@ref) entry is relative to its
    **current** mole fraction, since its activity is ``x_k\gamma_k``. For an
    end-member sitting at the solver's lower bound that is a statement about a
    vanishing phase, not about whether the solid solution would form.

## See also

  - [Thermochemistry](@ref sec-theory-thermo) — where ``\ln a`` enters ``\mu``
  - [Activity models](@ref sec-theory-activity) — the aqueous half
  - [Chemical Equilibrium](@ref sec-equilibrium) — declaring and solving with
    solid solutions
