# Thermodynamical functions

```@index
Pages = ["thermo_functions.md"]
```

```@docs
Callable
ThermoFunction
ThermoFunction(
    expr::Union{Symbol,Expr},
    params=Pair[],
    vars=[:T, :P, :t];
    ref=[:T => 298.15u"K", :P => 1u"bar", :t => 0u"s"],
)
thermo_function_library
∂
∫
calculate_molar_mass
apply
```
