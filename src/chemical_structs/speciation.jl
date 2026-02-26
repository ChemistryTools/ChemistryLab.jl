"""
    same_components(::AbstractVector{<:AbstractSpecies}) -> Function

Return the function to extract components from species vectors.

Returns `atoms_charge` for Species vectors, `oxides_charge` for CemSpecies vectors.
"""
same_components(::AbstractVector{<:AbstractSpecies}) = atoms_charge
same_components(::AbstractVector{<:CemSpecies}) = oxides_charge

"""
    item_order(::AbstractVector{<:AbstractSpecies}) -> Vector{Symbol}

Return the ordering vector for components.

Returns `ATOMIC_ORDER` for Species vectors, `OXIDE_ORDER` for CemSpecies vectors.
"""
item_order(::AbstractVector{<:AbstractSpecies}) = ATOMIC_ORDER
item_order(::AbstractVector{<:CemSpecies}) = OXIDE_ORDER

"""
    union_atoms(atom_dicts::AbstractVector{<:AbstractDict}, order_vec=ATOMIC_ORDER) -> Vector{Symbol}

Compute the union of all keys from dictionaries, sorted by a given order.

# Arguments

  - `atom_dicts`: vector of dictionaries (e.g., atomic compositions).
  - `order_vec`: ordering vector for sorting keys (default: ATOMIC_ORDER).

# Returns

  - Sorted vector of unique symbols appearing in any dictionary.

# Examples

```jldoctest
julia> d1 = OrderedDict(:H => 2, :O => 1);

julia> d2 = OrderedDict(:C => 1, :O => 2);

julia> union_atoms([d1, d2], ATOMIC_ORDER)
3-element Vector{Symbol}:
 :C
 :H
 :O
```
"""
function union_atoms(atom_vecs::AbstractVector{<:AbstractVector}, order_vec=ATOMIC_ORDER)
    function sortfunc(k)
        idx = findfirst(==(k), order_vec)
        return isnothing(idx) ? max(1, length(order_vec) - 1) : idx
    end
    return sort!(collect(union(atom_vecs...)); by=sortfunc)
end

union_atoms(atom_dicts::AbstractVector{<:AbstractDict}, order_vec=ATOMIC_ORDER) =
   union_atoms(collect.(keys.(atom_dicts)), order_vec)

union_atoms(species_list::AbstractVector{<:AbstractSpecies}, order_vec=ATOMIC_ORDER) =
   union_atoms(collect.(keys.(atoms.(species_list))), order_vec)

function idx_speciation(
    species_list, atoms_list::AbstractVector{Symbol};
    aggregate_state=[AS_AQUEOUS, AS_CRYSTAL, AS_GAS, AS_UNDEF],
    class=[SC_AQSOLUTE, SC_AQSOLVENT, SC_COMPONENT, SC_GASFLUID, SC_UNDEF],
    exclude_species=[],
    include_species=[],
)
    mask1 = [all(k .∈ Ref(atoms_list)) for k in keys.(atoms.(species_list))]
    mask2 = ChemistryLab.aggregate_state.(species_list) .∈ Ref(aggregate_state)
    mask3 = ChemistryLab.class.(species_list) .∈ Ref(class)
    mask4 = (species_list .∉ Ref(exclude_species)) .&&
            (symbol.(species_list) .∉ Ref(exclude_species))
    mask5 = (species_list .∈ Ref(include_species)) .||
            (symbol.(species_list) .∈ Ref(include_species))
    return (mask1 .&& mask2 .&& mask3 .&& mask4) .|| mask5
end

function speciation(
    species_list, atoms_list::AbstractVector{Symbol}; kwargs...
)
    return @view species_list[idx_speciation(species_list, atoms_list; kwargs...)]
end

function speciation(
    species_list, short_species_list::AbstractVector{<:AbstractSpecies}; kwargs...
)
    return speciation(species_list, union_atoms(short_species_list);
                 include_species = short_species_list, kwargs...)
end


function speciation(
    species_list, short_species_list_symbols::AbstractVector{<:AbstractString}; kwargs...
)
    short_species_list = @view species_list[symbol.(species_list) .∈ Ref(short_species_list_symbols)]
    return speciation(species_list, short_species_list; kwargs...)
end
