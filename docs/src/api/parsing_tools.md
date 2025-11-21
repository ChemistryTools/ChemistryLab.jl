# Parsing tools

```@index
Pages = ["parsing_tools.md"]
```

```@docs
stoich_coef_round
phreeqc_to_unicode
unicode_to_phreeqc
colored_formula
```

```@example colored_formula
using ChemistryLab # hide
for s in ("HSO₄²⁻", "HSO4-2", "(CaO)₄//₃(SiO₂)₁(H₂O)₁₃//₆")
    println(colored_formula(s))
end
```

```@docs
parse_formula
extract_charge
to_mendeleev
parse_equation
colored_equation
```

```@example colored_formula
println(colored_equation("Ca₄Al₂(OH)₁₄(H₂O)₆ + 6H⁺ = 2AlO₂⁻ + 16H₂O@ + 4Ca²⁺"))
```

```@docs
format_equation
```
