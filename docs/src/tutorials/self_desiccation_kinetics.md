# [Self-desiccation in time: the arrest as a result](@id sec-self-desiccation-kinetics)

[Self-desiccation](@ref sec-self-desiccation) takes Powers' 0.42 apart into a
water budget and reads the arrest off it at an assumed humidity. The degree of
hydration is an **input** there, and the page says so. This one integrates
instead: a rate law that reads the internal relative humidity of the paste it is
hydrating, so the arrest comes out of the trajectory.

What remains empirical is the **criterion** — [`humidity_factor`](@ref)
implements the cut that stops hydration below about 80 % RH — not the arrest
itself. That is the honest ceiling of a 0D framework, and
[the water budget](@ref sec-theory-water-budget) §4 explains why nothing in
thermodynamics can replace it.

## 1. The loop that closes

Three pieces already existed and were never connected:

| piece | what it gives |
|:--|:--|
| the rate law | consumes water, so the pore solution shrinks as the reaction runs |
| [`PoreHumidity`](@ref) | the saturation of the remaining pore space, and through a retention law its relative humidity |
| [`humidity_factor`](@ref) | the rate correction that vanishes as the humidity falls |

Composed, they are a feedback: hydration dries the paste, drying slows
hydration, and the reaction stops when the two balance. Nothing imposes where.

```@example sdk
using ChemistryLab
using OrdinaryDiffEq
using DynamicQuantities
using OrderedCollections
using Printf

substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
selection = split(
    "C3S C2S C3A C4AF Portlandite Jennite ettringite monosulphate12 C3AH6 C3FH6 H2O@"
)
cs = ChemicalSystem(
    speciation(substances, selection; aggregate_state = [AS_AQUEOUS]), CEMDATA_PRIMARIES
)
sp(name) = cs[name]

# Baroghel-Bouny et al. (1999), Table 5, mix CO — the same measured curve the
# static page uses, so the two are comparable.
law = VanGenuchten(; a = 37.5479e6, m = 1 / 2.1684)
nothing # hide
```

The clinker and the reactions are those of
[Cement clinker hydration kinetics](@ref); the only difference
is the `humidity` keyword.

```@example sdk
const COMPOSITION = (C3S = 0.619, C2S = 0.165, C3A = 0.08, C4AF = 0.087)

function build(wc; humidity)
    st = ChemicalState(cs)
    for (name, frac) in pairs(COMPOSITION)
        set_quantity!(st, string(name), frac * u"kg")
    end
    set_quantity!(st, "H2O@", wc * u"kg")

    # `PoreHumidity` needs the fresh paste as its volume reference: the pore
    # space it reports the saturation of is the one that existed at mixing.
    h = humidity === :pore ? PoreHumidity(law, cs; reference = st) : humidity

    specs = (
        ("C3S", PK84_PARAMS_C3S,
         OrderedDict(sp("C3S") => 1.0, sp("H2O@") => 103 / 30),
         OrderedDict(sp("Jennite") => 1.0, sp("Portlandite") => 4 / 3)),
        ("C2S", PK84_PARAMS_C2S,
         OrderedDict(sp("C2S") => 1.0, sp("H2O@") => 73 / 30),
         OrderedDict(sp("Jennite") => 1.0, sp("Portlandite") => 1 / 3)),
        ("C3A", PK84_PARAMS_C3A,
         OrderedDict(sp("C3A") => 1.0, sp("H2O@") => 6.0),
         OrderedDict(sp("C3AH6") => 1.0)),
        ("C4AF", PK84_PARAMS_C4AF,
         OrderedDict(sp("C4AF") => 1.0, sp("Portlandite") => 2.0, sp("H2O@") => 10.0),
         OrderedDict(sp("C3AH6") => 1.0, sp("C3FH6") => 1.0)),
    )
    rxns = AbstractReaction[]
    for (nm, pk, reac, prod) in specs
        rx = Reaction(reac, prod; symbol = nm)
        # α_max = 1.0: NO Powers cap. Whatever arrest appears is the humidity's.
        rx[:rate] = parrot_killoh_avrami(
            pk, nm; α_max = 1.0, blaine = 380.0u"m^2/kg", humidity = h
        )
        push!(rxns, rx)
    end
    kp = KineticsProblem(cs, rxns, st, (0.0, 90 * 86400.0))
    ks = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-10)
    return kp, integrate(kp, ks), h
end

α_end(sol, kp) = let d = degrees_of_hydration(sol, kp; times = [sol.t[end]])
    v = d["C3S"]
    v isa AbstractVector ? v[end] : v
end
nothing # hide
```

## 2. The arrest, and where it lands

```@example sdk
println("  w/c   α(C3S) coupled   α(C3S) uncoupled   Powers w/c/0.42   S at 90 d   RH at 90 d")
for wc in (0.25, 0.30, 0.35, 0.40, 0.50)
    kp, sol, h = build(wc; humidity = :pore)
    kpf, solf, _ = build(wc; humidity = nothing)
    n_end = ustrip.(us"mol", state_at(sol, kp, sol.t[end]).n)
    @printf("%6.2f  %14.4f  %18.4f  %17.4f  %11.4f  %10.4f\n",
            wc, α_end(sol, kp), α_end(solf, kpf), powers_alpha_max(wc),
            pore_saturation(h, n_end), h(n_end))
end
```

Three things to read off that table.

**The arrest is a result.** The uncoupled column is the same number at every
`w/c` — the rate law alone knows nothing about how much water there is, so it
reaches the same degree of hydration in ninety days whatever the mix. With the
humidity coupled, the degree of hydration depends on `w/c`, and does so because
the paste dried itself.

**It stops where `humidity_factor` cuts, at the saturation the static budget
assumes.** For the drier mixes the internal humidity lands on 0.800 and the pore
saturation on 0.786 — and 0.786 is exactly the ``S^\ast`` the
[static page](@ref sec-self-desiccation) reads off the *same* retention curve at
RH 0.80. The trajectory arrives at the number the budget assumes, by a different
route, with nothing arranged to make it so. That agreement is the closest thing
to a validation available here.

**The proportionality is gone, and that is informative.** The static
construction gives ``\alpha_{\max} \propto w/c`` exactly, which the static page
is careful to call structural and no evidence. Integrated, the ratio is not
constant: ``k = (w/c)/\alpha_{\max}`` comes out near 0.45 at the drier end and
rises with `w/c`, because a wet paste has not finished self-desiccating within
ninety days — the arrest is then set by the integration window rather than by
the water. Powers' 0.42 sits just below the dry-end value.

## 3. What was in the way

This run could not be performed before, and the reason is worth recording
because it is invisible from the object's own behavior. Every van Genuchten
retention law has an **unbounded** ``\mathrm{d}p_c/\mathrm{d}S`` at full
saturation, and a sealed paste starts saturated. An implicit solver
differentiates the rate law at its first step, so it received an infinite
Jacobian entry and returned at ``t = 0`` with every degree of hydration at zero.

[`PoreHumidity`](@ref) now returns the saturated value as a constant — water at
zero suction has ``a_w = 1``, which is not in doubt — so the derivative there is
zero instead of infinite. That is a **regularization and not an identity**: the
true derivative is unbounded. It applies within ``10^{-10}`` of saturation, in
practice the initial condition alone, and the solver leaves that point on its
first successful step.

## 4. What this does and does not predict

  - **predicted**: that a sealed paste arrests; that the arrest degree of
    hydration rises with `w/c`; the whole trajectory `α(t)`; and the pore
    saturation at which it stops, which agrees with the static budget;
  - **assumed**: the humidity threshold in [`humidity_factor`](@ref), which is
    empirical, and the retention curve, which is measured on one material at one
    age and applied at every degree of hydration;
  - **not represented**: transport. The reason a real paste stops is that the
    liquid path between its remaining water and an unreacted grain is broken,
    and a 0D model has no geometry — see
    [the water budget](@ref sec-theory-water-budget) §4. What this page does is
    put the empirical consequence of that fact inside a rate law, where it can
    at least act on the trajectory.

## 5. Reproducing this

Every block above runs when this page is built. See also
[Self-desiccation](@ref sec-self-desiccation) for the static budget,
[`PoreHumidity`](@ref), [`humidity_factor`](@ref),
[`parrot_killoh_avrami`](@ref) and [`VanGenuchten`](@ref).
