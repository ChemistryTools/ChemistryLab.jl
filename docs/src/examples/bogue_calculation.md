# Bogue Calculation

The way in which species and cementitious species are constructed in ChemistryLab and expressed as a linear combination of reference species opens the door to equilibrium calculations. It also makes it quite natural to retrieve Bogue's formulas and use them simply.

Bogue's formulas allow us to find the masses of $\ce{C3S}$, $\ce{C2S}$, $\ce{C3A}$ and $\ce{C4AF}$ as a function of the oxides ($\ce{CaO}$, $\ce{SiO2}$, $\ce{Al2O3}$ and $\ce{Fe2O3}$) that are regularly found in manufacturers' cement data sheets. However, using the `StoichMatrix` functions performs a molar decomposition of the species that we wish to decompose as a function of reference species. It is therefore possible to express the anhydrous of the cement as a function of the oxides in a cement data sheet.

## Via `ChemicalSystem`

The recommended approach builds a `ChemicalSystem` once, which computes both stoichiometric matrices and all derived structures automatically. The four clinker phases are the species; the four oxide components (C=$\ce{CaO}$, S=$\ce{SiO2}$, A=$\ce{Al2O3}$, F=$\ce{Fe2O3}$) are the primaries:

```@example Bogue
using ChemistryLab #hide
using DynamicQuantities #hide
using PrettyTables #hide
cemspecies = CemSpecies.(split("C3S C2S C3A C4AF"))
oxides = CemSpecies.(split("C S A F"))

cs = ChemicalSystem(cemspecies, oxides)
A = cs.SM.A
```

The stoichiometric matrix `cs.SM.A` (primaries × species) gives the molar decomposition of each clinker phase into oxides. Bogue's formulas are obtained by converting to mass and inverting:

```@example Bogue
# Molar mass of anhydrous phases
Mw  = map(x -> ustrip(us"g/mol", x.M), cemspecies)
# Molar mass of each oxide
Mwo = map(x -> ustrip(us"g/mol", x.M), oxides)
Aoa = Mwo .* A .* inv.(Mw)'

pprint(inv(Aoa), cemspecies, oxides; label=:name)
```

By taking a cement sheet with classic percentages of oxides ($\ce{CaO}$=65.6%, $\ce{SiO2}$=21.5%, $\ce{Al2O3}$=5.2% and $\ce{Fe2O3}$=2.8%), we obtain the anhydrous mass fractions of the cementitious material:

```@example Bogue
inv(Aoa) * [65.6, 21.5, 5.2, 2.8]
```

---

## Via `mass_matrix`

A more direct route uses `mass_matrix` on the canonical stoichiometric matrix. This expresses each species as a linear combination of the reference oxide components directly in mass rather than in moles, making Bogue's formulas immediate:

```@example Bogue
massCSM = mass_matrix(CanonicalStoichMatrix(cemspecies))
pprint(inv(massCSM.A), cemspecies, oxides; label=:name)
```
