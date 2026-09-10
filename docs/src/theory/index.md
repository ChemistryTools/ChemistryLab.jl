# [Theory](@id sec-theory)

This chapter answers *why this is the right calculation*, and its pages are
meant to be readable without running anything. The
[Tutorials](@ref sec-equilibrium) drive one feature at a time and the Examples
work a real case end to end; here the question is what the equations are, where
they come from, and where they stop being true.

## The three layers, and what each one assumes

ChemistryLab computes in three layers, and almost every surprise comes from a
layer's assumptions being carried into a regime it was not built for.

| layer | what it computes | what it assumes |
|:--|:--|:--|
| **chemical description** | formulas, species, reactions, the conservation matrix | nothing physical — this is bookkeeping, and it is exact |
| **equilibrium** | the composition minimizing the Gibbs energy under a conservation budget | one well-mixed phase per aggregate state, ideal molar volumes, and an **activity model** |
| **kinetics** | a trajectory in time, optionally re-equilibrating the solution at every step | a rate law per reaction, and that the rate law's arguments are available |

The activity model is where the second layer stops being ideal, so it is the
first thing to read and the first thing to suspect:
[Activity models](@ref sec-theory-activity).

## Reading order

Start with [Activity models](@ref sec-theory-activity), which every equilibrium
calculation depends on whether or not it is mentioned. The derivation of the
minimization itself, its optimality certificate and the constraint machinery are
in [Chemical Equilibrium](@ref sec-equilibrium) for now; the cement-specific
theory — the water budget, Powers' coefficient, and what a 0D calculation can
and cannot predict about a drying paste — is in
[Self-desiccation](@ref sec-self-desiccation).
