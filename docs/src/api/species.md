# Species

```@index
Pages = ["species.md"]
```

```@docs
AbstractSpecies
AggregateState
Class
Species
Base.isequal(s1::AbstractSpecies, s2::AbstractSpecies)
charge(s::AbstractSpecies)
CemSpecies
name(s::AbstractSpecies)
symbol(s::AbstractSpecies)
formula(s::AbstractSpecies)
cemformula
atoms(s::AbstractSpecies)
atoms_charge(s::AbstractSpecies)
oxides
oxides_charge
components(s::Species)
components(s::CemSpecies)
aggregate_state(s::AbstractSpecies)
class(s::AbstractSpecies)
with_class
properties(s::AbstractSpecies)
check_mendeleev(s::AbstractSpecies)
apply(func::Function, s::S, args...; kwargs...) where {S<:AbstractSpecies}
```
