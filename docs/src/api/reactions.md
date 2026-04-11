# Reactions

```@index
Pages = ["reactions.md"]
```

```@docs
Reaction
CemReaction
symbol(r::Reaction)
equation(r::Reaction)
colored(r::Reaction)
equal_sign(r::Reaction)
reactants
products
charge
properties(r::Reaction)
simplify_reaction
apply(func::Function, r::Reaction, args...; kwargs...)
```
