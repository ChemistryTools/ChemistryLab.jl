# Bogue Calculation

The way in which species and cementitious species are constructed in ChemistryLab and expressed as a linear combination of reference species opens the door to equilibrium calculations. It also makes it quite natural to retrieve Bogue's formulas and use them simply.

Bogue's formulas allow us to find the masses of $\text{C}_3\text{S}$, $\text{C}_2\text{S}$, $\text{C}_3\text{A}$ and $\text{C}_4\text{AF}$ as a function of the oxides ($\text{CaO}$, $\text{SiO}_2$, $\text{Al}_2\text{O}_3$ and $\text{Fe}_2\text{O}_3$) that are regularly found in manufacturers' cement data sheets. However, using the `StoichMatrix` functions performs a molar decomposition of the species that we wish to decompose as a function of reference species. It is therefore possible to express the anhydrous of the cement as a function of the oxides in a cement data sheet.

```@example Bogue
using ChemistryLab #hide
using DynamicQuantities #hide
cemspecies = CemSpecies.(split("C3S C2S C3A C4AF"))
oxides = CemSpecies.(split("C S A F"))
SM = StoichMatrix(cemspecies, oxides)
A = SM.A

using PrettyTables #hide
```

Bogue's formulas are thus found by converting species into mass and inverting the matrix.

```@example Bogue
# Molar mass of anhdrous phases
Mw = map(x -> ustrip(us"g/mol", x.M), cemspecies)
# Molar mass of each oxide
Mwo = map(x -> ustrip(us"g/mol", x.M), oxides)
Aoa = Mwo .* A .* inv.(Mw)'

pprint(inv(Aoa), cemspecies, oxides; label=:name)
```

By taking a cement sheet with a classic percentages of oxides ($\text{CaO}$=65.6%, $\text{SiO}_2$=21.5%, $\text{Al}_2\text{O}_3$=5.2% and $\text{Fe}_2\text{O}_3$=2.8%), we then obtain the anhydrous masses of the cementitious material. 

```@example Bogue
inv(Aoa) *  [65.6, 21.5, 5.2, 2.8]
```

---

Another, faster way to do this is to take advantage of the `mass=true` option of the `CanonicalStoichMatrix` and `StoichMatrix` functions. This option allows species to be expressed as a linear combination of the reference species in mass rather than in moles.

```@example Bogue
massCSM = mass_matrix(CanonicalStoichMatrix(cemspecies))
```

Bogue's formulas are then immediate.

```@example Bogue
pprint(inv(massCSM.A), cemspecies, oxides; label=:name)
```
