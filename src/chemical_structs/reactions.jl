"""
    fwd_arrows

Collection of forward arrow symbols used in chemical reaction notation.

Contains symbols representing forward reaction directions, including:

  - Simple arrows: '>', '→'
  - Various Unicode arrows with different styles and weights
  - Specialized arrows for reaction notation

Used to define reaction directionality from reactants to products.
"""
const fwd_arrows = ['>', '→', '↣', '↦', '⇾', '⟶', '⟼', '⥟', '⥟', '⇀', '⇁', '⇒', '⟾']

"""
    bwd_arrows

Collection of backward arrow symbols used in chemical reaction notation.

Contains symbols representing reverse reaction directions, including:

  - Simple arrows: '<', '←'
  - Various Unicode arrows with different styles and weights
  - Specialized arrows for reverse reaction notation

Used to define reaction directionality from products to reactants.
"""
const bwd_arrows = ['<', '←', '↢', '↤', '⇽', '⟵', '⟻', '⥚', '⥞', '↼', '↽', '⇐', '⟽']

"""
    double_arrows

Collection of double arrow symbols representing equilibrium reactions.

Contains symbols representing:

  - Simple equilibrium: '↔'
  - Various Unicode double arrows
  - Specialized equilibrium symbols
  - Bidirectional reaction indicators

Used to denote reversible reactions and equilibrium states.
"""
const double_arrows = ['↔', '⟷', '⇄', '⇆', '⇌', '⇋', '⇔', '⟺']

"""
    pure_rate_arrows

Collection of specialized arrows for rate-based reaction notation.

Contains symbols commonly used to represent:

  - Reaction rates
  - Kinetic directions
  - Specialized reaction mechanisms

These are often used in more advanced chemical kinetics notation.
"""
const pure_rate_arrows = ['⇐', '⟽', '⇒', '⟾', '⇔', '⟺']

"""
    equal_signs

Collection of equality signs used in chemical reaction equations.

Contains various forms of equality operators including:

  - Standard equals sign: '='
  - Definition operators: '≔'
  - Specialized equality symbols
  - Assignment operators

Used to separate reactants from products in balanced equations.
"""
const equal_signs = ['=', '≔', '⩴', '≕']

const EQUAL_REACTION = vcat(
    fwd_arrows, bwd_arrows, double_arrows, pure_rate_arrows, equal_signs
)
const EQUAL_REACTION_SET = Set(EQUAL_REACTION)

"""
    struct Reaction{SR<:AbstractSpecies,TR<:Number,SP<:AbstractSpecies,TP<:Number}

Representation of a chemical reaction with reactants and products.

# Fields

  - `equation::String`: Unicode equation string.
  - `colored::String`: colored terminal representation.
  - `reactants::OrderedDict{SR,TR}`: species => coefficient for reactants.
  - `products::OrderedDict{SP,TP}`: species => coefficient for products.
  - `equal_sign::Char`: equality operator character.
  - `properties::OrderedDict{Symbol,PropertyType}`: thermodynamic and other properties.

# Examples

```julia
julia> length(reactants(r))
r = Reaction("2H2 + O2 = 2H2O");

julia> length(products(r))
2
```
"""
struct Reaction{SR<:AbstractSpecies,TR<:Number,SP<:AbstractSpecies,TP<:Number}
    equation::String
    colored::String
    reactants::OrderedDict{SR,TR}
    products::OrderedDict{SP,TP}
    equal_sign::Char
    properties::OrderedDict{Symbol,PropertyType}
end

"""
    equation(r::Reaction) -> String

Return the equation string of the reaction.

# Examples

```julia
julia> equation(r)
r = Reaction("H2 + O2 = H2O");
```
"""
equation(r::Reaction) = r.equation

"""
    colored(r::Reaction) -> String

Return the colored terminal representation of the reaction.

# Examples

```julia
julia> colored(r)  # Returns string with ANSI color codes
r = Reaction("H2 + O2 = H2O");
```
"""
colored(r::Reaction) = r.colored

"""
    reactants(r::Reaction) -> OrderedDict

Return the reactants dictionary (species => coefficient).

# Examples

```julia
julia> reactants(r)
r = Reaction("2H2 + O2 = 2H2O");
```
"""
reactants(r::Reaction) = r.reactants

"""
    products(r::Reaction) -> OrderedDict

Return the products dictionary (species => coefficient).

# Examples

```julia
julia> products(r)
r = Reaction("2H2 + O2 = 2H2O");
```
"""
products(r::Reaction) = r.products

"""
    equal_sign(r::Reaction) -> Char

Return the equality operator character of the reaction.

# Examples

```julia
julia> equal_sign(r)
r = Reaction("H2 + O2 = H2O");
```
"""
equal_sign(r::Reaction) = r.equal_sign

"""
    properties(r::Reaction) -> OrderedDict{Symbol,PropertyType}

Return the properties dictionary of the reaction.

# Examples

```julia
julia> properties(r)
r = Reaction("H2 + O2 = H2O");
```
"""
properties(r::Reaction) = r.properties

"""
    Base.getindex(r::Reaction, i::Symbol) -> Any

Access a reaction property by symbol key.
Return `nothing` if the property is not found.

# Examples

```julia
julia> r[:ΔrH⁰]
r = Reaction("H2 + O2 = H2O");

julia> r[:nonexistent]
r[:ΔrH⁰] = -241.8;
```
"""
Base.getindex(r::Reaction, i::Symbol) = get(properties(r), i, nothing)

"""
    Base.getindex(r::Reaction, s::AbstractSpecies) -> Number

Get the stoichiometric coefficient for a species in the reaction.
Return negative values for reactants, positive for products, and 0 if the species is not present.

# Examples

```julia
julia> r[h2]
r = Reaction("2H2 + O2 = 2H2O");

julia> r[o2]
h2 = Species("H2");

julia> r[h2o]
-2

julia> r[co2]
o2 = Species("O2");
```
"""
function Base.getindex(r::Reaction, s::AbstractSpecies)
    coef = get(r.products, s, nothing)
    if isnothing(coef)
        coef = get(r.reactants, s, nothing)
        if !isnothing(coef)
            coef = -coef
        else
            return 0
        end
    end
    return coef
end

"""
    Base.setindex!(r::Reaction, value, i::Symbol)

Set a property value for the reaction.

# Examples

```julia
julia> r[:ΔrH⁰]
r = Reaction("H2 + O2 = H2O");
```
"""
Base.setindex!(r::Reaction, value, i::Symbol) = setindex!(properties(r), value, i)

"""
    Base.getproperty(r::Reaction, sym::Symbol) -> Any

Access reaction fields or registered properties.
Throws an error if the symbol is neither a field nor a property.

# Examples

```julia
julia> r.ΔrH⁰
r = Reaction("H2 + O2 = H2O");

julia> r.equation
r[:ΔrH⁰] = -241.8;
```
"""
function Base.getproperty(r::Reaction, sym::Symbol)
    if sym in fieldnames(typeof(r))
        return getfield(r, sym)
    elseif sym in keys(properties(r))
        return properties(r)[sym]
    else
        error("Symbol '$sym' is neither a field nor a registered property.")
    end
end

"""
    Base.haskey(r::Reaction, sym::Symbol) -> Bool

Check if a property key exists in the reaction properties dictionary.

# Examples

```julia
julia> haskey(r, :ΔrH⁰)
r = Reaction("H2 + O2 = H2O");

julia> haskey(r, :nonexistent)
r[:ΔrH⁰] = -241.8;
```
"""
function Base.haskey(r::Reaction, sym::Symbol)
    return haskey(properties(r), sym)
end

"""
    Base.setproperty!(r::Reaction, sym::Symbol, value)

Set a property value, preventing direct modification of structural fields.

# Examples

```julia
julia> r.ΔrH⁰
r = Reaction("H2 + O2 = H2O");
```
"""
function Base.setproperty!(r::Reaction, sym::Symbol, value)
    if !ismissing(value)
        if sym in fieldnames(typeof(r))
            error(
                "Cannot modify field '$sym' directly. Use constructor or dedicated methods."
            )
        else
            properties(r)[sym] = value
        end
    end
    return r
end

"""
    Base.iterate(r::Reaction, state=(1, nothing))

Iterate over all species in the reaction with signed coefficients.
Yields (species, coefficient) pairs where coefficients are negative for reactants
and positive for products.

# Examples

```julia
julia> collect(r)
r = Reaction("2H2 + O2 = 2H2O");
```
"""
function Base.iterate(r::Reaction, state=(1, nothing))
    idx, inner_state = state
    if idx == 1
        if inner_state === nothing
            inner_state = iterate(reactants(r))
        else
            inner_state = iterate(reactants(r), inner_state)
        end
        if inner_state === nothing
            return iterate(r, (2, nothing))
        else
            (k, v), new_state = inner_state
            return (k, -v), (1, new_state)
        end
    elseif idx == 2
        if inner_state === nothing
            inner_state = iterate(products(r))
        else
            inner_state = iterate(products(r), inner_state)
        end
        if inner_state === nothing
            return nothing
        else
            kv, new_state = inner_state
            return kv, (2, new_state)
        end
    else
        return nothing
    end
end

function Base.length(r::Reaction)
    return length(reactants(r)) + length(products(r))
end

"""
    Base.keys(r::Reaction)

Return an iterator over all species in the reaction (reactants and products).

# Examples

```julia
julia> collect(keys(r))
r = Reaction("2H2 + O2 = 2H2O");
```
"""
function Base.keys(r::Reaction)
    return Iterators.flatten((keys(reactants(r)), keys(products(r))))
end

"""
    Base.values(r::Reaction)

Return an iterator over all stoichiometric coefficients (negative for reactants, positive for products).

# Examples

```julia
julia> collect(values(r))
r = Reaction("2H2 + O2 = 2H2O");
```
"""
function Base.values(r::Reaction)
    vals1 = (-v for v in values(reactants(r)))
    vals2 = values(products(r))
    return Iterators.flatten((vals1, vals2))
end

"""
    remove_zeros(d::AbstractDict) -> AbstractDict

Remove all entries with zero values from a dictionary.

# Examples

```julia
julia> remove_zeros(d)
d = OrderedDict("H2" => 2, "O2" => 0, "H2O" => 1);
```
"""
function remove_zeros(d::AbstractDict)
    for (k, v) in d
        if iszero(v)
            delete!(d, k)
        end
    end
    return d
end

"""
    complete_thermo_functions(r::Reaction)

Compute reaction thermodynamic properties from species properties.
Calculates ΔrCp⁰, ΔrS⁰, ΔrH⁰, ΔrG⁰, and ΔrV if all species have the required properties.

# Examples

```julia
julia> r.ΔrCp⁰
h2 = Species("H2"); h2.Cp⁰ = 28.8;
```
"""
function complete_thermo_functions(r::Reaction)
    species_list = keys(r)
    if !isempty(species_list)
        if all(x -> haskey(x, :Cp⁰), species_list)
            r.ΔrCp⁰ = sum(ν * s.Cp⁰ for (s, ν) in r)
        end
        if all(x -> haskey(x, :S⁰), species_list)
            r.ΔrS⁰ = sum(ν * s.S⁰ for (s, ν) in r)
        end
        if all(x -> haskey(x, :ΔfH⁰), species_list)
            r.ΔrH⁰ = sum(ν * s.ΔfH⁰ for (s, ν) in r)
        end
        if all(x -> haskey(x, :ΔfG⁰), species_list)
            r.ΔrG⁰ = sum(ν * s.ΔfG⁰ for (s, ν) in r)
        end
        if all(x -> haskey(x, :Vm), species_list)
            r.ΔrV = sum(ν * s.Vm for (s, ν) in r)
        end
        r.charge = sum(ν * charge(s) for (s, ν) in r)
    else
        r.charge = 0
    end
end

"""
    Reaction(equation::AbstractString, S::Type{<:AbstractSpecies}=Species; properties, side, species_list) -> Reaction

Construct a Reaction from an equation string.

# Arguments

  - `equation`: reaction equation string (e.g., "2H2 + O2 = 2H2O").
  - `S`: species type to use (default: Species).
  - `properties`: property dictionary (default: empty OrderedDict).
  - `side`: how to split species - :none, :sign, :reactants, :products (default: :none).
  - `species_list`: optional list of known species for lookup.

# Examples

```julia
julia> length(reactants(r))
r = Reaction("H2 + 0.5O2 = H2O");

julia> length(products(r))
2
```
"""
function Reaction(
    equation::AbstractString,
    S::Type{<:AbstractSpecies}=Species;
    properties::AbstractDict=OrderedDict{Symbol,PropertyType}(),
    side::Symbol=:none,
    species_list=nothing,
)
    reactants, products, equal_sign = parse_equation(equation)
    r = Reaction(
        equation,
        colored_equation(equation),
        ordered_dict_with_default(
            (
                find_species(k, species_list) => stoich_coef_round(v) for
                (k, v) in reactants if !iszero(v) && !startswith(k, "Zz") && !startswith(k, "e")
            ),
            S,
            Number,
        ),
        ordered_dict_with_default(
            (
                find_species(k, species_list) => stoich_coef_round(v) for
                (k, v) in products if !iszero(v) && !startswith(k, "Zz") && !startswith(k, "e")
            ),
            S,
            Number,
        ),
        equal_sign,
        OrderedDict{Symbol,PropertyType}(properties),
    )
    complete_thermo_functions(r)
    if side == :none
        return r
    else
        return Reaction(r; side=side)
    end
end

"""
    CemReaction(equation::AbstractString, args...; kwargs...) -> Reaction

Construct a Reaction using CemSpecies from an equation string.
Convenience constructor equivalent to `Reaction(equation, CemSpecies, args...; kwargs...)`.

# Examples

```julia
julia> length(reactants(r))
r = CemReaction("CaO + H2O = Ca(OH)2");

julia> length(products(r))
2
```
"""
function CemReaction(equation::AbstractString, args...; kwargs...)
    Reaction(equation, CemSpecies, args...; kwargs...)
end

"""
    split_species_by_stoich(species_stoich::AbstractDict{S,T}; side=:sign) where {S<:AbstractSpecies,T<:Number} -> (OrderedDict, OrderedDict)

Split a species-coefficient dictionary into reactants and products.

# Arguments

  - `species_stoich`: dictionary mapping species to signed stoichiometric coefficients.
  - `side`: splitting criterion - :sign (by coefficient sign), :reactants, :left, :products, :right.

# Returns

  - Tuple of (reactants_dict, products_dict) with positive coefficients.

# Examples

```julia
julia> split_species_by_stoich(s)
h2 = Species("H2"); o2 = Species("O2"); h2o = Species("H2O");
```
"""
function split_species_by_stoich(
    species_stoich::AbstractDict{S,T}; side::Symbol=:sign
) where {S<:AbstractSpecies,T<:Number}
    reactants = OrderedDict{S,T}()
    products = OrderedDict{S,T}()
    for (species, coef) in species_stoich
        if !iszero(coef)
            if try
                side == :reactants || side == :left || (coef < 0 && side == :sign)
            catch
                false
            end
                reactants[species] = -stoich_coef_round(coef)
            else
                products[species] = stoich_coef_round(coef)
            end
        end
    end
    return reactants, products
end

"""
    merge_species_by_stoich(reactants::AbstractDict{SR,TR}, products::AbstractDict{SP,TP}) where {SR,TR,SP,TP} -> OrderedDict

Merge reactants and products into a single dictionary with signed coefficients.
Reactants get negative coefficients, products get positive coefficients.

# Examples

```julia
julia> merge_species_by_stoich(r, p)
r = OrderedDict(Species("H2") => 2);
```
"""
function merge_species_by_stoich(
    reactants::AbstractDict{SR,TR}, products::AbstractDict{SP,TP}
) where {SR<:AbstractSpecies,TR<:Number,SP<:AbstractSpecies,TP<:Number}
    return merge(
        +,
        ordered_dict_with_default(
            (species => -stoich_coef_round(coef) for (species, coef) in reactants), SR, TR
        ),
        ordered_dict_with_default(
            (species => stoich_coef_round(coef) for (species, coef) in products), SP, TP
        ),
    )
end

"""
    format_side(side::AbstractDict{S,T}) where {S<:AbstractSpecies,T<:Number} -> (String, String, Int)

Format one side of a reaction equation.

# Arguments

  - `side`: dictionary of species => coefficient for one side of the reaction.

# Returns

  - Tuple of (equation_string, colored_string, total_charge).

# Examples

```julia
julia> format_side(s)
h2 = Species("H2"); o2 = Species("O2");
```
"""
function format_side(side::AbstractDict{S,T}) where {S<:AbstractSpecies,T<:Number}
    equation = String[]
    coleq = String[]
    ch = 0
    Zz = root_type(S)("Zz")
    for (species, coef) in side
        if !iszero(coef) && species != Zz
            coeff_str = isone(coef) ? "" : string(stoich_coef_round(coef))
            coeff_str = add_parentheses_if_needed(coeff_str)
            coeff_str = replace(coeff_str, " " => "", "*" => "")
            push!(equation, coeff_str * unicode(species))
            push!(coleq, string(COL_STOICH_EXT(coeff_str)) * colored(species))
            ch += coef * charge(species)
        end
    end
    if isempty(equation)
        equation = "∅"
        coleq = "∅"
    end
    return join(equation, " + "), join(coleq, " + "), ch
end

"""
    Reaction(reactants::AbstractDict{SR,TR}, products::AbstractDict{SP,TP}; equal_sign='=', properties, side) where {SR,TR,SP,TP} -> Reaction

Construct a Reaction from reactants and products dictionaries.

# Arguments

  - `reactants`: dictionary mapping reactant species to coefficients.
  - `products`: dictionary mapping product species to coefficients.
  - `equal_sign`: equality operator character (default '=').
  - `properties`: property dictionary (default: empty OrderedDict).
  - `side`: how to reorganize species - :none, :sign, :reactants, :products (default: :none).
    Automatically balances electron charges in the equation.

# Examples

```julia
julia> length(products(r))
h2 = Species("H2"); o2 = Species("O2"); h2o = Species("H2O");
```
"""
function Reaction(
    reactants::AbstractDict{SR,TR},
    products::AbstractDict{SP,TP};
    equal_sign='=',
    properties::AbstractDict=OrderedDict{Symbol,PropertyType}(),
    side::Symbol=:none,
) where {SR<:AbstractSpecies,TR<:Number,SP<:AbstractSpecies,TP<:Number}
    if side ∈ (:sign, :products, :right, :reactants, :left)
        reactants, products = split_species_by_stoich(
            merge_species_by_stoich(reactants, products); side=side
        )
    end
    delete!(reactants, root_type(SR)("Zz"))
    delete!(reactants, root_type(SR)("e"))
    delete!(products, root_type(SP)("Zz"))
    delete!(products, root_type(SP)("e"))
    sreac, creac, charge_left = format_side(reactants)
    sprod, cprod, charge_right = format_side(products)
    charge_diff = charge_right - charge_left
    if !isapprox(charge_diff, 0; atol=1e-4)
        needed_e = if charge_diff < 0
            -stoich_coef_round(charge_diff)
        else
            stoich_coef_round(charge_diff)
        end
        e_term = needed_e == 1 ? "e⁻" : "$needed_e" * "e⁻"
        ce_term = if needed_e == 1
            "e⁻"
        else
            string(COL_STOICH_EXT(add_parentheses_if_needed("$needed_e"))) * "e⁻"
        end
        if charge_diff < 0
            sreac = isempty(sreac) ? e_term : "$sreac + $e_term"
            creac = isempty(creac) ? e_term : "$creac + $ce_term"
        else
            sprod = isempty(sprod) ? e_term : "$sprod + $e_term"
            cprod = isempty(cprod) ? e_term : "$cprod + $ce_term"
        end
    end
    equation = sreac * " " * string(equal_sign) * " " * sprod
    colored = creac * " " * string(COL_PAR(string(equal_sign))) * " " * cprod
    r = Reaction(
        equation,
        colored,
        OrderedDict{SR,TR}(reactants),
        OrderedDict{SP,TP}(products),
        equal_sign,
        OrderedDict{Symbol,PropertyType}(properties),
    )
    complete_thermo_functions(r)
    return r
end

"""
    Reaction(species_stoich::AbstractDict{S,T}; equal_sign='=', properties, side=:sign) where {S,T} -> Reaction

Construct a Reaction from a dictionary with signed stoichiometric coefficients.

# Arguments

  - `species_stoich`: dictionary mapping species to signed coefficients (negative = reactants, positive = products).
  - `equal_sign`: equality operator character (default '=').
  - `properties`: property dictionary (default: empty OrderedDict).
  - `side`: splitting criterion (default: :sign).

# Examples

```julia
julia> length(reactants(r))
h2 = Species("H2"); o2 = Species("O2"); h2o = Species("H2O");
```
"""
function Reaction(
    species_stoich::AbstractDict{S,T};
    equal_sign::Char='=',
    properties::AbstractDict=OrderedDict{Symbol,PropertyType}(),
    side::Symbol=:sign,
) where {S<:AbstractSpecies,T<:Number}
    reactants, products = split_species_by_stoich(species_stoich; side=side)
    return Reaction(
        reactants,
        products;
        equal_sign=equal_sign,
        properties=OrderedDict{Symbol,PropertyType}(properties),
    )
end

"""
    Base.convert(::Type{Reaction}, s::S) where {S<:AbstractSpecies} -> Reaction

Convert a species to a trivial Reaction (species = species).

# Examples

```julia
julia> r[h2o]
h2o = Species("H2O");
```
"""
function Base.convert(::Type{Reaction}, s::S) where {S<:AbstractSpecies}
    Reaction(OrderedDict(s => 1))
end

"""
    Base.convert(::Type{Reaction{U,T}}, s::S) where {U,T,S} -> Reaction

Convert a species to a typed Reaction.

# Examples

```julia
julia> r[h2o]
h2o = Species("H2O");
```
"""
function Base.convert(
    ::Type{Reaction{U,T}}, s::S
) where {U<:AbstractSpecies,T<:Number,S<:AbstractSpecies}
    Reaction(OrderedDict(s => 1))
end

"""
    Reaction(s::S) where {S<:AbstractSpecies} -> Reaction

Construct a trivial Reaction from a single species.

# Examples

```julia
julia> r[h2o]
h2o = Species("H2O");
```
"""
Reaction(s::S) where {S<:AbstractSpecies} = Reaction(OrderedDict(s => 1))

"""
    Reaction{U,T}(s::S) where {U,T,S} -> Reaction

Construct a typed Reaction from a single species.

# Examples

```julia
julia> r[h2o]
h2o = Species("H2O");
```
"""
function Reaction{U,T}(s::S) where {U<:AbstractSpecies,T<:Number,S<:AbstractSpecies}
    Reaction(OrderedDict(s => 1))
end

"""
    Reaction(r::R; equal_sign, properties, side) where {R<:Reaction} -> Reaction

Copy constructor for Reaction with optional field overrides.

# Arguments

  - `r`: source Reaction.
  - `equal_sign`: override equality operator (default: keep original).
  - `properties`: override properties (default: keep original).
  - `side`: reorganization criterion (default: :none).

# Examples

```julia
julia> r2.equal_sign
r = Reaction("H2 + O2 = H2O");
```
"""
function Reaction(
    r::R; equal_sign=r.equal_sign, properties=r.properties, side::Symbol=:none
) where {R<:Reaction}
    Reaction(
        reactants(r),
        products(r);
        equal_sign=equal_sign,
        side=side,
        properties=OrderedDict{Symbol,PropertyType}(properties),
    )
end

"""
    simplify_reaction(r::Reaction) -> Reaction

Simplify a reaction by canceling common species from both sides.

# Examples

```julia
julia> length(reactants(rs))
h2o = Species("H2O");

julia> length(products(rs))
r = Reaction(OrderedDict(h2o => 2), OrderedDict(h2o => 1));
```
"""
function simplify_reaction(r::Reaction)
    reac = remove_zeros(deepcopy(reactants(r)))
    prod = remove_zeros(deepcopy(products(r)))
    common_species = intersect(keys(reac), keys(prod))
    for species in common_species
        coef = prod[species] - reac[species]
        if iszero(coef)
            delete!(reac, species)
            delete!(prod, species)
        elseif try
            coef > 0
        catch
            true
        end
            prod[species] = coef
            delete!(reac, species)
        else
            reac[species] = -coef
            delete!(prod, species)
        end
    end
    return Reaction(reac, prod; equal_sign=equal_sign(r), properties=properties(r))
end

"""
    scale_stoich!(species_stoich::AbstractDict{<:AbstractSpecies,<:Number})

Scale stoichiometric coefficients by their GCD if all are integers or rationals.
Modifies the dictionary in place to ensure integer coefficients when possible.

# Arguments

  - `species_stoich`: dictionary mapping species to stoichiometric coefficients

# Examples

```julia
julia> d
d = OrderedDict(Species("H2") => 2, Species("O2") => 1);

julia> d2
scale_stoich!(d);
```
"""
function scale_stoich!(species_stoich::AbstractDict{<:AbstractSpecies,<:Number})
    v = values(species_stoich)
    if all(x -> x isa Integer || x isa Rational, v)
        mult = gcd([numerator(x) for x in v]...)
        for k in keys(species_stoich)
            species_stoich[k] *= mult
        end
    end
end

"""
    build_species_stoich(species::AbstractVector{<:AbstractSpecies}; scaling=1, auto_scale=false) -> OrderedDict

Build stoichiometric coefficients from a species vector using stoichiometric matrix analysis.
The first species is treated as the dependent component.

# Arguments

  - `species`: vector of species (first is the dependent component)
  - `scaling`: scaling factor for all coefficients (default: 1)
  - `auto_scale`: if true, scale by GCD to get integer coefficients (default: false)

# Returns

  - OrderedDict mapping species to signed stoichiometric coefficients (negative for reactants)

# Examples

```julia
julia> build_species_stoich([h2o, h2, o2])
h2 = Species("H2"); o2 = Species("O2"); h2o = Species("H2O");

julia> build_species_stoich([h2o, h2, o2]; auto_scale=true)
OrderedDict(h2o => -1, h2 => 1, o2 => 0.5)
```
"""
function build_species_stoich(
    species::AbstractVector{<:AbstractSpecies}; scaling=1, auto_scale=false
)
    A, indep_comp, dep_comp = stoich_matrix(
        species[1:1], species[2:end]; display=false, involve_all_atoms=true
    )
    S, T = promote_type(typeof.(indep_comp)..., typeof.(dep_comp)...), eltype(A)
    species_stoich = OrderedDict{S,T}()
    species_stoich[dep_comp[1]] = -scaling
    for (i, s) in enumerate(indep_comp)
        species_stoich[s] = A[i, 1] * scaling
    end
    if auto_scale
        scale_stoich!(species_stoich)
    end
    return species_stoich
end

"""
    Reaction(species::AbstractVector{<:AbstractSpecies}; equal_sign='=', properties, scaling=1, auto_scale=false, side=:sign) -> Reaction

Construct a balanced Reaction from a vector of species.
The first species is treated as the dependent component, and stoichiometric
coefficients are computed automatically.

# Arguments

  - `species`: vector of species to balance (first is dependent component)
  - `equal_sign`: equality operator character (default: '=')
  - `properties`: property dictionary (default: empty OrderedDict)
  - `scaling`: scaling factor for all coefficients (default: 1)
  - `auto_scale`: if true, scale by GCD (default: false)
  - `side`: splitting criterion (default: :sign)

# Returns

  - A balanced Reaction object

# Examples

```julia
julia> r[h2]
h2 = Species("H2"); o2 = Species("O2"); h2o = Species("H2O");

julia> r[o2]
r = Reaction([h2o, h2, o2]);

julia> r2[h2]
1

julia> r2[o2]
0.5
```
"""
function Reaction(
    species::AbstractVector{<:AbstractSpecies};
    equal_sign='=',
    properties::AbstractDict=OrderedDict{Symbol,PropertyType}(),
    scaling=1,
    auto_scale=false,
    side::Symbol=:sign,
)
    species_stoich = build_species_stoich(species; scaling=scaling, auto_scale=auto_scale)
    return Reaction(
        species_stoich;
        equal_sign=equal_sign,
        properties=OrderedDict{Symbol,PropertyType}(properties),
        side=side,
    )
end

"""
    Reaction(reac::AbstractVector{<:AbstractSpecies}, prod::AbstractVector{<:AbstractSpecies}; kwargs...) -> Reaction

Construct a balanced Reaction from separate reactant and product vectors.
Stoichiometric coefficients are computed automatically to balance the reaction.

# Arguments

  - `reac`: vector of reactant species
  - `prod`: vector of product species
  - `equal_sign`: equality operator character (default: '=')
  - `properties`: property dictionary (default: empty OrderedDict)
  - `scaling`: scaling factor for all coefficients (default: 1)
  - `auto_scale`: if true, scale by GCD (default: false)
  - `side`: splitting criterion (default: :none)

# Returns

  - A balanced Reaction object

# Examples

```julia
julia> r[h2]
h2 = Species("H2"); o2 = Species("O2"); h2o = Species("H2O");

julia> r[o2]
r = Reaction([h2, o2], [h2o]);

julia> r[h2o]
-1

julia> r2[h2]
-0.5

julia> r2[o2]
1

julia> r2[h2o]
r2 = Reaction([h2, o2], [h2o]; auto_scale=true);
```
"""
function Reaction(
    reac::AbstractVector{<:AbstractSpecies},
    prod::AbstractVector{<:AbstractSpecies};
    equal_sign='=',
    properties::AbstractDict=OrderedDict{Symbol,PropertyType}(),
    scaling=1,
    auto_scale=false,
    side::Symbol=:none,
)
    species = [reac; prod]
    species_stoich = build_species_stoich(species; scaling=scaling, auto_scale=auto_scale)
    S, T = keytype(species_stoich), valtype(species_stoich)
    if side != :none
        return Reaction(
            species_stoich;
            equal_sign=equal_sign,
            properties=OrderedDict{Symbol,PropertyType}(properties),
            side=side,
        )
    else
        return Reaction(
            ordered_dict_with_default(
                (k => -v for (k, v) in species_stoich if k in reac), S, T
            ),
            ordered_dict_with_default(
                (k => v for (k, v) in species_stoich if k in prod), S, T
            );
            equal_sign=equal_sign,
            properties=OrderedDict{Symbol,PropertyType}(properties),
            side=:none,
        )
    end
end

"""
    *(ν::Number, s::AbstractSpecies) -> Reaction

Create a Reaction with a single species and stoichiometric coefficient.

# Arguments

  - `ν`: stoichiometric coefficient
  - `s`: species to include in the reaction

# Returns

  - A Reaction object with the single species and given coefficient

# Examples

```julia
julia> r[h2o]
h2o = Species("H2O");

julia> r2[h2o]
r = 2 * h2o;
```
"""
*(ν::Number, s::AbstractSpecies) = Reaction(OrderedDict(s => ν))

"""
    *(ν::Number, r::Reaction) -> Reaction

Multiply all stoichiometric coefficients in a reaction by a scalar.

# Arguments

  - `ν`: scaling factor
  - `r`: reaction to scale

# Returns

  - A new Reaction with all coefficients multiplied by ν

# Examples

```julia
julia> r2[o2]
r = Reaction("H2 + 0.5O2 = H2O");

julia> r2[h2]
r2 = 2 * r;
```
"""
function *(
    ν::Number, r::Reaction{SR,TR,SP,TP}
) where {SR<:AbstractSpecies,TR<:Number,SP<:AbstractSpecies,TP<:Number}
    Reaction(
        ordered_dict_with_default((k => ν * v for (k, v) in reactants(r)), SR, TR),
        ordered_dict_with_default((k => ν * v for (k, v) in products(r)), SP, TP);
        equal_sign=r.equal_sign,
        properties=r.properties,
    )
end

"""
    -(s::AbstractSpecies) -> Reaction

Create a Reaction with a single species with coefficient -1.

# Arguments

  - `s`: species to include in the reaction

# Returns

  - A Reaction object with the single species and coefficient -1

# Examples

```julia
julia> r[h2o]
h2o = Species("H2O");
```
"""
-(s::AbstractSpecies) = Reaction(OrderedDict(s => -1))

"""
    -(r::Reaction) -> Reaction

Reverse a reaction (swap reactants and products).

# Arguments

  - `r`: reaction to reverse

# Returns

  - A new Reaction with reactants and products swapped

# Examples

```julia
julia> rr[h2]
r = Reaction("H2O = H2 + 0.5O2");

julia> rr[o2]
rr = -r;
```
"""
-(r::Reaction) = Reaction(
    products(r), reactants(r); equal_sign=r.equal_sign, properties=r.properties
)

"""
    +(s::S1, t::S2) where {S1<:AbstractSpecies,S2<:AbstractSpecies} -> Reaction

Add two species to create a Reaction.

# Arguments

  - `s`: first species
  - `t`: second species

# Returns

  - A Reaction with both species as reactants (coefficient 1 each)

# Examples

```julia
julia> r[h2]
h2 = Species("H2"); o2 = Species("O2");

julia> r[o2]
r = h2 + o2;

julia> r2[h2]
1
```
"""
function +(s::S1, t::S2) where {S1<:AbstractSpecies,S2<:AbstractSpecies}
    S = promote_type(S1, S2)
    s == t ? Reaction(OrderedDict(S(s) => 2)) : Reaction(OrderedDict(S(s) => 1, S(t) => 1))
end

"""
    -(s::S1, t::S2) where {S1<:AbstractSpecies,S2<:AbstractSpecies} -> Reaction

Subtract two species to create a Reaction.

# Arguments

  - `s`: first species (positive coefficient)
  - `t`: second species (negative coefficient)

# Returns

  - A Reaction with s as reactant and t as product

# Examples

```julia
julia> r[h2]
h2 = Species("H2"); o2 = Species("O2");

julia> r[o2]
r = h2 - o2;

julia> length(reactants(r2))
1
```
"""
function -(s::S1, t::S2) where {S1<:AbstractSpecies,S2<:AbstractSpecies}
    S = promote_type(S1, S2)
    if s == t
        Reaction(OrderedDict{S,Number}())
    else
        Reaction(OrderedDict(S(s) => 1, S(t) => -1))
    end
end

"""
    add_stoich(d1::AbstractDict{S1,T1}, d2::AbstractDict{S2,T2}) where {S1<:AbstractSpecies,T1<:Number,S2<:AbstractSpecies,T2<:Number} -> OrderedDict

Add stoichiometric coefficients from two dictionaries.

# Arguments

  - `d1`: first dictionary of species => coefficients
  - `d2`: second dictionary of species => coefficients

# Returns

  - A new dictionary with combined coefficients

# Examples

```julia
julia> add_stoich(d1, d2)
d1 = OrderedDict(Species("H2") => 2);

julia> add_stoich(d1, d3)
d2 = OrderedDict(Species("O2") => 1);
```
"""
function add_stoich(
    d1::AbstractDict{S1,T1}, d2::AbstractDict{S2,T2}
) where {S1<:AbstractSpecies,T1<:Number,S2<:AbstractSpecies,T2<:Number}
    S = promote_type(S1, S2)
    T = promote_type(T1, T2)
    d = OrderedDict{S,T}()
    for (k, v) in d1
        d[k] = get(d, k, 0) + v
    end
    for (k, v) in d2
        d[k] = get(d, k, 0) + v
    end
    return d
end

"""
    +(r::R, s::S) where {R<:Reaction,S<:AbstractSpecies} -> Reaction

Add a species to a reaction.

# Arguments

  - `r`: reaction to modify
  - `s`: species to add as product

# Returns

  - A new Reaction with the species added as product

# Examples

```julia
julia> r2[o2]
r = Reaction("H2 = H2O");

julia> length(products(r2))
o2 = Species("O2");
```
"""
function +(r::R, s::S) where {R<:Reaction,S<:AbstractSpecies}
    return Reaction(
        reactants(r),
        add_stoich(products(r), OrderedDict(s => 1));
        equal_sign=r.equal_sign,
        properties=r.properties,
    )
end

"""
    -(r::R, s::S) where {R<:Reaction,S<:AbstractSpecies} -> Reaction

Subtract a species from a reaction.

# Arguments

  - `r`: reaction to modify
  - `s`: species to remove from products (or add as reactant)

# Returns

  - A new Reaction with the species subtracted

# Examples

```julia
julia> r2[h2]
r = Reaction("H2 + O2 = H2O");
```
"""
function -(r::R, s::S) where {R<:Reaction,S<:AbstractSpecies}
    return Reaction(
        reactants(r),
        add_stoich(products(r), OrderedDict(s => -1));
        equal_sign=r.equal_sign,
        properties=r.properties,
    )
end

+(s::S, r::R) where {S<:AbstractSpecies,R<:Reaction} = +(r, s)
-(s::S, r::R) where {S<:AbstractSpecies,R<:Reaction} = +(s, -r)

"""
    +(r::R, u::U) where {R<:Reaction,U<:Reaction} -> Reaction

Add two reactions.

# Arguments

  - `r`: first reaction
  - `u`: second reaction

# Returns

  - A new Reaction combining both reactions

# Examples

```julia
julia> length(reactants(r3))
r1 = Reaction("H2 = 2H");

julia> length(products(r3))
r2 = Reaction("O2 = 2O");
```
"""
function +(r::R, u::U) where {R<:Reaction,U<:Reaction}
    return Reaction(
        add_stoich(reactants(r), reactants(u)),
        add_stoich(products(r), products(u));
        equal_sign=r.equal_sign,
        properties=merge(properties(r), properties(u)),
    )
end

"""
    -(r::R, u::U) where {R<:Reaction,U<:Reaction} -> Reaction

Subtract two reactions.

# Arguments

  - `r`: first reaction
  - `u`: second reaction to subtract

# Returns

  - A new Reaction representing r - u

# Examples

```julia
julia> r3[Species("H2")]
r1 = Reaction("2H2 + O2 = 2H2O");

julia> r3[Species("H2O")]
r2 = Reaction("H2 = H2O");
```
"""
function -(r::R, u::U) where {R<:Reaction,U<:Reaction}
    return Reaction(
        add_stoich(reactants(r), products(u)),
        add_stoich(products(r), reactants(u));
        equal_sign=r.equal_sign,
        properties=merge(properties(r), properties(u)),
    )
end

"""
    EQUAL_OPS

Union of all supported equality operators for chemical reactions.

Combines forward arrows, backward arrows, double arrows, rate arrows,
and equality signs (excluding the first element of each collection to avoid duplicates).

This collection is used to dynamically generate reaction operator methods.

  - Forward: →, ↣, ↦, ⇾, ⟶, ⟼, ⥟, ⇀, ⇁, ⇒, ⟾
  - Backward: ←, ↢, ↤, ⇽, ⟵, ⟻, ⥚, ⥞, ↼, ↽, ⇐, ⟽
  - Equilibrium: ↔, ⟷, ⇄, ⇆, ⇌, ⇋, ⇔, ⟺
  - Rate: ⇐, ⟽, ⇒, ⟾
  - Equality: ≔, ⩴, ≕
"""
const EQUAL_OPS = union(
    fwd_arrows[2:end],
    bwd_arrows[2:end],
    double_arrows,
    pure_rate_arrows,
    equal_signs[2:end],
)

for OP in Symbol.(EQUAL_OPS)
    @eval begin
        $OP(r, s) = Reaction(
            -Reaction(r) + Reaction(s); equal_sign=first(string($OP)), side=:sign
        )
    end

    @eval @doc """
        $($OP)(r::Reaction, s::Reaction)

    Dynamically generated reaction operator method for chemical equations using the '$($OP)' symbol.

    This method:
    1. Takes two Reaction objects (r and s)
    2. Reverses the first reaction (r → -r)
    3. Adds the second reaction (s)
    4. Creates a new Reaction with:
        - The combined species from both reactions
        - The equality operator set to '$($OP)'
        - Species split by sign (reactants vs products)

    This allows writing chemical equations with natural notation.

    # Arguments
    - `r`: First reaction (will be reversed in the operation)
    - `s`: Second reaction (will be added as-is)

    # Returns
    - A new Reaction object combining both reactions with the '$($OP)' operator

    # Examples
    ```julia
    julia> h2 = Species("H2"); o2 = Species("O2"); h2o = Species("H2O");
    julia> r1 = Reaction([h2, o2], [h2o]);  # H2 + O2 → H2O
    julia> r2 = Reaction([h2o], [h2, o2]);  # H2O → H2 + O2

    julia> result = $($OP)(r1, r2);  # Creates combined reaction with '$($OP)' operator
    julia> equal_sign(result)
    '$($OP)'
    ```

    # See also
    All supported operators: →, ↣, ↦, ⇾, ⟶, ⟼, ⥟, ⇀, ⇁, ⇒, ⟾, ←, ↢, ↤, ⇽, ⟵, ⟻, ⥚, ⥞, ↼, ↽, ⇐, ⟽, ↔, ⟷,     ⇄, ⇆, ⇌, ⇋, ⇔, ⟺, ≔, ⩴, ≕
    """ $OP
end

"""
    Base.show(io::IO, r::Reaction)

Display a reaction in a compact form.

# Arguments

  - `io`: output stream
  - `r`: reaction to display

# Examples

```julia
julia> print(r)
r = Reaction("H2 + O2 = H2O");
```
"""
function Base.show(io::IO, r::Reaction)
    print(io, colored(r))
end

"""
    Base.show(io::IO, ::MIME"text/plain", r::Reaction)

Display a reaction in a detailed form.

# Arguments

  - `io`: output stream
  - `r`: reaction to display

# Examples

```julia
julia> show(stdout, MIME"text/plain"(), r)
r = Reaction("H2 + O2 = H2O");
```
"""
function Base.show(io::IO, ::MIME"text/plain", r::Reaction)
    println(io, colored(r))
    pad = 10
    if length(reactants(r)) > 0
        println(
            io,
            lpad("reactants", pad),
            ": ",
            join(["$(colored(k)) => $v" for (k, v) in reactants(r)], ", "),
        )
    else
        println(io, lpad("reactants", pad), ": ∅")
    end
    pr = length(properties(r)) > 0 ? println : print
    if length(products(r)) > 0
        pr(
            io,
            lpad("products", pad),
            ": ",
            join(["$(colored(k)) => $v" for (k, v) in products(r)], ", "),
        )
    else
        pr(io, lpad("products", pad), ": ∅")
    end
    if length(properties(r)) > 0
        print(
            io,
            lpad("properties", pad),
            ": ",
            join(["$k = $v" for (k, v) in properties(r)], "\n" * repeat(" ", pad + 2)),
        )
    end
end

"""
    apply(func::Function, r::Reaction{SR,TR,SP,TP}, args...; kwargs...) where {SR<:AbstractSpecies,TR<:Number,SP<:AbstractSpecies,TP<:Number}

Apply a function to all species and coefficients in a reaction.

# Arguments

  - `func`: function to apply to species and coefficients
  - `r`: reaction to transform
  - `args...`: additional arguments for func
  - `kwargs...`: additional keyword arguments

# Returns

  - A new Reaction with transformed species and coefficients

# Examples

```julia
julia> r_prime.equation
r = Reaction("H2 + O2 = H2O");

julia> collect(keys(r_prime))
r_prime = apply(s -> uppercase(name(s)), r);
```
"""
function apply(
    func::Function, r::Reaction{SR,TR,SP,TP}, args...; kwargs...
) where {SR<:AbstractSpecies,TR<:Number,SP<:AbstractSpecies,TP<:Number}
    tryfunc(v) =
        if v isa Quantity
            (
                try
                    func(ustrip(v), args...; kwargs...) *
                    func(dimension(v), args...; kwargs...)
                catch
                    try
                        func(ustrip(v), args...; kwargs...) * dimension(v)
                    catch
                        v
                    end
                end
            )
        else
            (
                try
                    func(v, args...; kwargs...)
                catch
                    v
                end
            )
        end
    reac = OrderedDict{SR,TR}(
        apply(func, k, args...; kwargs...) => tryfunc(v) for (k, v) in reactants(r)
    )
    prod = OrderedDict{SP,TP}(
        apply(func, k, args...; kwargs...) => tryfunc(v) for (k, v) in products(r)
    )
    newReaction = Reaction(
        reac,
        prod;
        equal_sign=get(kwargs, :equal_sign, equal_sign(r)),
        properties=OrderedDict{Symbol,PropertyType}(
            k => v for (k, v) in get(kwargs, :properties, properties(r))
        ),
        side=get(kwargs, :side, :none),
    )
    for (k, v) in properties(r)
        newReaction[k] = tryfunc(v)
    end
    return newReaction
end
