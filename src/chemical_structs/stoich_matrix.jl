"""
    same_components(::Vector{<:AbstractSpecies}) -> Function

Return the function to extract components from species vectors.

Returns `atoms_charge` for Species vectors, `oxides_charge` for CemSpecies vectors.
"""
same_components(::Vector{<:AbstractSpecies}) = atoms_charge
same_components(::Vector{<:CemSpecies}) = oxides_charge

"""
    item_order(::Vector{<:AbstractSpecies}) -> Vector{Symbol}

Return the ordering vector for components.

Returns `ATOMIC_ORDER` for Species vectors, `OXIDE_ORDER` for CemSpecies vectors.
"""
item_order(::Vector{<:AbstractSpecies}) = ATOMIC_ORDER
item_order(::Vector{<:CemSpecies}) = OXIDE_ORDER

"""
    union_atoms(atom_dicts::Vector{<:AbstractDict}, order_vec=ATOMIC_ORDER) -> Vector{Symbol}

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
function union_atoms(atom_dicts::Vector{<:AbstractDict}, order_vec=ATOMIC_ORDER)
    function sortfunc(k)
        idx = findfirst(==(k), order_vec)
        return isnothing(idx) ? max(1, length(order_vec)-1) : idx
    end
    sort!(collect(union(keys.(atom_dicts)...)); by=sortfunc)
end

"""
    pprint(A::AbstractMatrix, indep_comp_names::Vector, dep_comp_names::Vector)

Print a stoichiometric matrix with colored formatting.

# Arguments

  - `A`: stoichiometric matrix.
  - `indep_comp_names`: row labels (independent components).
  - `dep_comp_names`: column labels (dependent components).

Uses text highlighters to color positive (red), negative (blue), and zero (concealed) values.
"""
function pprint(
    A::AbstractMatrix, indep_comp_names::Vector, dep_comp_names::Vector
)
    hl_p = TextHighlighter((data, i, j) -> (data[i, j] > 0), crayon"bold light_red")
    hl_n = TextHighlighter((data, i, j) -> (data[i, j] < 0), crayon"bold light_blue")
    hl_z = TextHighlighter((data, i, j) -> (data[i, j] == 0), crayon"conceal")
    try
        pretty_table(
            A;
            column_labels=dep_comp_names,
            row_labels=indep_comp_names,
            highlighters=[hl_p, hl_n, hl_z],
            style=TextTableStyle(;
                row_label=crayon"magenta bold",
                first_line_column_label=crayon"cyan bold",
                table_border=crayon"green bold",
            ),
        )
    catch
        pretty_table(
            A;
            column_labels=dep_comp_names,
            row_labels=indep_comp_names,
            style=TextTableStyle(;
                row_label=crayon"magenta bold",
                first_line_column_label=crayon"cyan bold",
                table_border=crayon"green bold",
            ),
        )
    end
end

"""
    stoich_matrix_to_equations(A::AbstractMatrix, indep_comp_names::AbstractVector, dep_comp_names::AbstractVector; scaling=1, pprint=true, equal_sign='=') -> Vector{String}

Convert a stoichiometric matrix to chemical equation strings.

# Arguments

  - `A`: stoichiometric matrix.
  - `indep_comp_names`: independent component names.
  - `dep_comp_names`: dependent component names.
  - `scaling`: scaling factor for coefficients (default 1).
  - `pprint`: if true, print equations to stdout (default true).
  - `equal_sign`: equality operator character (default '=').

# Returns

  - Vector of equation strings.

Each column of A represents a dependent component expressed as a linear combination
of independent components.
"""
function stoich_matrix_to_equations(
    A::AbstractMatrix,
    indep_comp_names::AbstractVector,
    dep_comp_names::AbstractVector;
    scaling=1,
    pprint=false,
    equal_sign='=',
)
    eqns = String[]
    pad = 11
    for (j, sp) in enumerate(dep_comp_names)
        if sp in indep_comp_names
            if pprint
                println(
                    rpad("$(sp)", pad),
                    "| $(colored_formula(sp)) $(string(COL_PAR(string(equal_sign)))) $(colored_formula(sp))",
                )
            end
        else
            coeffs = OrderedDict(zip(indep_comp_names, -A[:, j]))
            coeffs[sp] = 1
            eqn = format_equation(coeffs; scaling=scaling, equal_sign=equal_sign)
            push!(eqns, eqn)
            if pprint
                println(rpad("$(sp)", pad), "| ", colored_equation(eqn))
            end
        end
    end
    return eqns
end

"""
    stoich_matrix_to_reactions(A::AbstractMatrix, indep_comp::AbstractVector{<:AbstractSpecies}, dep_comp::AbstractVector{<:AbstractSpecies}; scaling=1, pprint=true, equal_sign='=') -> Vector{Reaction}

Convert a stoichiometric matrix to Reaction objects.

# Arguments

  - `A`: stoichiometric matrix.
  - `indep_comp`: independent component species.
  - `dep_comp`: dependent component species.
  - `scaling`: scaling factor for coefficients (default 1).
  - `pprint`: if true, print reactions to stdout (default true).
  - `equal_sign`: equality operator character (default '=').

# Returns

  - Vector of Reaction objects.
"""
function stoich_matrix_to_reactions(
    A::AbstractMatrix,
    indep_comp::AbstractVector{<:AbstractSpecies},
    dep_comp::AbstractVector{<:AbstractSpecies};
    scaling=1,
    pprint=true,
    equal_sign='=',
)
    eqns = Reaction[]
    pad = 11
    for (j, sp) in enumerate(dep_comp)
        if sp in indep_comp
            if pprint
                println(
                    rpad("$(symbol(sp))", pad),
                    "│ $(colored(sp)) $(string(COL_PAR(string(equal_sign)))) $(colored(sp))",
                )
            end
        else
            coeffs = OrderedDict(zip([sp; indep_comp], [1; -A[:, j]]))
            eqn = scaling * Reaction(coeffs)
            push!(eqns, eqn)
            if pprint
                println(rpad("$(symbol(sp))", pad), "│ ", colored(eqn))
            end
        end
    end
    return eqns
end

"""
    canonical_stoich_matrix(species::Vector{<:AbstractSpecies}; pprint=false, label=:symbol, mass=false) -> (Matrix, Vector{Symbol})

Build the canonical stoichiometric matrix for a species list.

# Arguments

  - `species`: vector of species to analyze.
  - `pprint`: if true, print the matrix (default true).
  - `label`: field name for labeling columns (default :symbol).
  - `mass`: if true, compute mass-based matrix (default false).

# Returns

  - Tuple of (stoichiometric_matrix, involved_atoms_vector).

The matrix has rows for each atom/oxide and columns for each species.
Entry (i,j) is the coefficient of atom i in species j.

# Examples

```jldoctest
julia> species = [Species("H2O"), Species("H2"), Species("O2")];

julia> A, atoms = canonical_stoich_matrix(species);

julia> size(A)
(2, 3)
```
"""
function canonical_stoich_matrix(
    species::Vector{<:AbstractSpecies}; pprint=false, label=:symbol, mass=false
)
    involved_atoms_dicts = same_components(species).(species)
    involved_atoms = union_atoms(involved_atoms_dicts, item_order(species))
    T = promote_type(valtype.(involved_atoms_dicts)...)
    A = zeros(T, length(involved_atoms), length(species))
    for (j, atoms) in enumerate(involved_atoms_dicts)
        for (i, atom) in enumerate(involved_atoms)
            A[i, j] = get(atoms, atom, zero(T))
        end
    end

    if mass
        S = promote_type(root_type.(typeof.(species))...)
        Msp = ustrip.(getproperty.(species, :M))
        Mat = ustrip.(getproperty.(S.(string.(involved_atoms)), :M))
        A = Mat .* A .* inv.(Msp)'
    end

    if pprint
        ChemistryLab.pprint(A, involved_atoms, eval(label).(species))
    end
    return A, involved_atoms
end

"""
    stoich_matrix(vs::Vector{<:AbstractSpecies}, candidate_primaries::Vector{<:AbstractSpecies}=vs; pprint=false, label=:symbol, involve_all_atoms=true, reorder_primaries=false, mass=false) -> (Matrix, Vector{AbstractSpecies}, Vector{AbstractSpecies})

Compute the stoichiometric matrix expressing dependent species in terms of independent components.

# Arguments

  - `vs`: vector of dependent species to express.
  - `candidate_primaries`: vector of candidate primary species (default: vs).
  - `pprint`: if true, print the matrix (default true).
  - `label`: field name for labeling (default :symbol).
  - `involve_all_atoms`: if true, include all atoms from candidates (default true).
  - `reorder_primaries`: if true, use QR pivoting to select primaries (default false).
  - `mass`: if true, compute mass-based matrix (default false).

# Returns

  - Tuple of (A, independent_components, dependent_components) where A[i,j] is the
    coefficient of independent component i in the decomposition of dependent component j.

# Examples

```jldoctest
julia> h2o = Species("H2O");
       h2 = Species("H2");
       o2 = Species("O2");

julia> A, indep, dep = stoich_matrix([h2o], [h2, o2]);

julia> size(A)
(2, 1)
```
"""
function stoich_matrix(
    vs::Vector{<:AbstractSpecies},
    candidate_primaries::Vector{<:AbstractSpecies}=vs;
    pprint=false,
    label=:symbol,
    involve_all_atoms=true,
    reorder_primaries=false,
    mass=false,
)
    safe_rank(A; rtol=1e-6) =
        try
            rank(A, rtol=rtol)
        catch
            try
                rank(A)
            catch
                min(size(A)...)
            end
        end
    safe_pinv(A) =
        try
            pinv(A)
        catch
            inv(A)
        end

    all_species = union(vs, candidate_primaries)
    vec_components = same_components(all_species)

    S = promote_type(typeof.(vs)..., typeof.(candidate_primaries)...)

    species = S[]
    append!(species, vs)
    num_initial_species = length(species)
    initial_involved_atoms = if involve_all_atoms
        union_atoms(vec_components.(all_species), item_order(species))
    else
        union_atoms(vec_components.(species), item_order(species))
    end
    candidate_primaries = deepcopy(candidate_primaries)

    for x in candidate_primaries
        idx = findfirst(y -> x == y, species)
        if isnothing(idx) &&
            all(k -> first(k) ∈ initial_involved_atoms || first(k) == :Zz, vec_components(x))
            push!(species, x)
        end
    end

    initial_involved_atoms = union_atoms(vec_components.(species), item_order(species))

    SpType(::Vector) = Species
    SpType(::Vector{<:CemSpecies}) = CemSpecies
    Zz = SpType(species)("Zz")
    charged = :Zz ∈ initial_involved_atoms
    if charged
        if Zz ∉ species
            push!(species, Zz)
        end
        if Zz ∉ candidate_primaries
            push!(candidate_primaries, Zz)
        end
    end

    M, involved_atoms = canonical_stoich_matrix(species; pprint=false)
    redox =
        charged && safe_rank(M[:, 1:(end - 1)]) != safe_rank(M[1:(end - 1), 1:(end - 1)])

    if !redox && charged
        pop!(species)
        M = M[1:(end - 1), 1:(end - 1)]
    end

    cols_candidates = [findfirst(y -> y == x, species) for x in candidate_primaries]
    filter!(x -> !isnothing(x), cols_candidates)
    M_subset = M[:, cols_candidates]
    # if size(M_subset, 1) >= size(M_subset, 2)
    #     independent_cols_indices = cols_candidates
    # else
    #     F = qr(M_subset, reorder_primaries ? Val(true) : NoPivot())
    #     r = Int(safe_rank(M_subset))
    #     pivot_idx = reorder_primaries ? F.p[1:r] : 1:r
    #     independent_cols_indices = sort(cols_candidates[pivot_idx])
    # end

    r = Int(safe_rank(M_subset))
    if reorder_primaries
        F = qr(M_subset, Val(true))
        pivot_idx = F.p[1:r]
        independent_cols_indices = sort(cols_candidates[pivot_idx])
    else
        pivot_idx = Int[]
        current_matrix = zeros(eltype(M_subset), size(M_subset, 1), 0)
        for j in axes(M_subset, 2)
            candidate = hcat(current_matrix, M_subset[:, j])
            if Int(safe_rank(candidate)) > size(current_matrix, 2)
                push!(pivot_idx, j)
                current_matrix = candidate
                if length(pivot_idx) == r
                    break
                end
            end
        end
        independent_cols_indices = cols_candidates[pivot_idx]
    end

    sort!(
        independent_cols_indices;
        by=x ->
            symbol(species[x]) !== "H2O@" &&
            symbol(species[x]) !== "H2O" &&
            symbol(species[x]) !== "H₂O" &&
            symbol(species[x]) !== "H",
    )
    M_indep = M[:, independent_cols_indices]
    M_indep = promote_type(typeof.(M_indep)...).(M_indep)
    A = stoich_coef_round.(safe_pinv(M_indep) * M)

    indep_comp = species[independent_cols_indices]
    dep_comp = species[1:num_initial_species]
    A = A[:, 1:num_initial_species]

    if redox && Zz ∈ dep_comp
        A = A[:, 1:(end - 1)]
        dep_comp = dep_comp[1:(end - 1)]
    end

    if mass
        S = promote_type(root_type.(typeof.(species))...)
        Mdep = ustrip.(getproperty.(dep_comp, :M))
        Mindep = ustrip.(getproperty.(indep_comp, :M))
        A = Mindep .* A .* inv.(Mdep)'
    end

    zero_rows = all(iszero.(A), dims=2)[:, 1]
    A = A[.!zero_rows, :]
    indep_comp = indep_comp[.!zero_rows]

    if pprint
        ChemistryLab.pprint(A, eval(label).(indep_comp), eval(label).(dep_comp))
    end

    return A, indep_comp, dep_comp
end

"""
    const oxides_as_species

Vector of Species representing cement oxides (C, S, A, F, etc.) in atomic composition.
"""
const oxides_as_species = [Species(d; symbol=string(k)) for (k, d) in CEMENT_TO_MENDELEEV]

"""
    const Aoxides

Canonical stoichiometric matrix for cement oxides.
"""
const Aoxides, atoms_in_oxides = canonical_stoich_matrix(oxides_as_species; pprint=false)

"""
    const order_atom_in_oxides

Dictionary mapping atoms to their row indices in the oxide stoichiometric matrix.
"""
const order_atom_in_oxides = Dict(atom => i for (i, atom) in enumerate(atoms_in_oxides))
