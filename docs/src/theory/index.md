# [Theory](@id sec-theory)

This chapter answers *why this is the right calculation*, and its pages are
meant to be readable without running anything. The
[Tutorials](@ref sec-equilibrium) drive one feature at a time and the Examples
work a real case end to end; here the question is what the equations are, where
they come from, and where they stop being true.

!!! note "What belongs on a theory page, and what does not"
    These pages carry the general theory: definitions, derivations, and the
    meaning of each term. Code appears only where the correspondence between a
    formula and its implementation *is* the point — a signature, a one-line
    expression. **A theory page does not build a chemical system, does not
    solve, and does not print a table of numbers**; that is what the
    [Applications](@ref sec-app-activity-models) are for, and every quantitative
    claim made here is measured there or asserted in the test suite. The
    division exists so that there is one place to look for each kind of thing.

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

**[Thermochemistry](@ref sec-theory-thermo)** first: it fixes the notation the
rest of the chapter uses, and that notation is the code's. Then the two places a
mixture stops being ideal — **[Activity models](@ref sec-theory-activity)** for
the aqueous phase, which every equilibrium calculation depends on whether or not
it is mentioned, and **[Solid solutions](@ref sec-theory-solid-solutions)** for a
phase of variable composition.

**[The water budget of a hydrating paste](@ref sec-theory-water-budget)** is the
cement-specific chapter, and the one to read if the question is why a
calculation predicts a threshold at ``w/c \approx 0.30`` where Powers reports
0.42. It is also where the limits of a 0D framework are argued rather than
asserted.

The derivation of the minimization itself, its optimality certificate and the
constraint machinery are still in
[Chemical Equilibrium](@ref sec-equilibrium) — theory living in a tutorial, and
a move waiting to happen.
