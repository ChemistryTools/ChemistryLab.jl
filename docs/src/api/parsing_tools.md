# Parsing tools

```@index
Pages = ["parsing_tools.md"]
```

## Examples

```@example colored_formula
using ChemistryLab # hide
for s in ("HSO₄²⁻", "HSO4-2", "(CaO)₄//₃(SiO₂)₁(H₂O)₁₃//₆")
    println(colored_formula(s))
end
```

```@example colored_formula
println(colored_equation("Ca₄Al₂(OH)₁₄(H₂O)₆ + 6H⁺ = 2AlO₂⁻ + 16H₂O@ + 4Ca²⁺"))
```

## Reference

```@autodocs
Modules = [ChemistryLab]
Pages = ["chemical_structs/parsing_tools.jl"]
```
