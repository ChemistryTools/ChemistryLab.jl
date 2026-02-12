# Formulas

```@index
Pages = ["formulas.md"]
```

```@docs
AtomGroup{T<:Number}
AtomGroup(sym::Symbol)
Base.convert(::Type{AtomGroup}, sym::Symbol)
stoichtype
Formula
expr(f::Formula)
phreeqc(f::Formula)
unicode(f::Formula)
colored(f::Formula)
composition(f::Formula)
charge(f::Formula)
check_mendeleev(f::Formula)
calculate_molar_mass
apply
pprint
```
