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
        StoichMatrix{T,P}

Container holding a stoichiometric matrix `A` together with the
`primaries` (independent components) and the full `species` vector.

# Fields

    - `A`: matrix (components × species) of stoichiometric coefficients.
    - `primaries`: vector of independent components (Symbols or `Species`).
    - `species`: vector of `Species` corresponding to the columns of `A`.
    - `N`: nullspace matrix (if possibly, built by the inner constructor): columns are vectors of the kernel of `A`.

# Examples

```julia
julia> species = [Species("H₂O"), Species("OH⁻"), Species("H⁺")]
3-element Vector{Species{Int64}}:
 H₂O {H₂O} [H₂O ◆ H2O]
 OH⁻ {OH⁻} [OH⁻ ◆ OH-]
 H⁺ {H⁺} [H⁺ ◆ H+]

julia> SM = StoichMatrix(species)
┌─────┬─────┬─────┬────┐
│     │ H₂O │ OH⁻ │ H⁺ │
├─────┼─────┼─────┼────┤
│ H₂O │   1 │     │  1 │
│ OH⁻ │     │   1 │ -1 │
└─────┴─────┴─────┴────┘
```
"""
struct StoichMatrix{T<:Number,P,V<:AbstractVector{P},M<:AbstractMatrix{T},S<:AbstractSpecies} <: AbstractMatrix{T}
    A::M
    primaries::V
    species::Vector{S}
    N::M
end

Base.eltype(::StoichMatrix{T}) where {T} = T

primtype(::StoichMatrix{T,P}) where {T,P} = P

function StoichMatrix(A::M, primaries::Union{Vector{Symbol},Vector{S}}, species::Vector{S}) where {M<:AbstractMatrix, S<:AbstractSpecies}
    indices_in = [findfirst(x -> x == p, species) for p in primaries]
    if any(x->isnothing(x), indices_in)
        return StoichMatrix(A, primaries, species, similar(A, 0, 0))
    else
        p, n = size(A)
        N = similar(A, n, n-p)
        indices_out = setdiff(eachindex(species), indices_in)
        N[indices_out, :] .= Diagonal(ones(eltype(A), length(indices_out)))
        N[indices_in, :] .= -A[:, indices_out]
        return StoichMatrix(A, primaries, species, N)
    end
end

for f in(:size, :length, :getindex, :ndims, :iterate)
    @eval Base.$f(SM::StoichMatrix,  args...; kwargs...) = $f(SM.A,  args...; kwargs...)
end

Base.show(io::IO, SM::StoichMatrix) = show(io, SM.A)

function Base.show(::IO, ::MIME"text/plain", SM::StoichMatrix)
    (; A, primaries, species) = SM
    column_labels = try symbol.(species) catch; species end
    row_labels = try symbol.(primaries) catch; primaries end
    formatters = [(v, i, j) -> iszero(v) ? "" : v]
    pretty_table(
        A;
        column_labels=column_labels,
        row_labels=row_labels,
        formatters=formatters,
    )
    return nothing
end

apply(f::Function, SM::StoichMatrix,  args...; kwargs...) =
    StoichMatrix(f(SM.A, args...; kwargs...), SM.primaries, SM.species)

"""
    pprint(A::AbstractMatrix, indep_comp_names::Vector, dep_comp_names::Vector)

Print a stoichiometric matrix with colored formatting.

# Arguments

  - `A`: stoichiometric matrix.
  - `indep_comp_names`: row labels (independent components).
  - `dep_comp_names`: column labels (dependent components).

Uses text highlighters to color positive (red), negative (blue), and zero (concealed) values.
"""
function pprint(A::AbstractMatrix, indep_comp_names::Vector, dep_comp_names::Vector; row_label = :symbol, col_label = :symbol, label = :identity)
    if label != :identity
        row_label = label
        col_label = label
    end
    column_labels = try eval(col_label).(dep_comp_names) catch; dep_comp_names end
    row_labels = try eval(col_label).(indep_comp_names) catch; indep_comp_names end
    hl_p = TextHighlighter((data, i, j) -> (data[i, j] > 0), crayon"bold light_red")
    hl_n = TextHighlighter((data, i, j) -> (data[i, j] < 0), crayon"bold light_blue")
    hl_z = TextHighlighter((data, i, j) -> (iszero(data[i, j])), crayon"conceal")
    try
        pretty_table(
            A;
            column_labels=column_labels,
            row_labels=row_labels,
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
            column_labels=column_labels,
            row_labels=row_labels,
            style=TextTableStyle(;
                row_label=crayon"magenta bold",
                first_line_column_label=crayon"cyan bold",
                table_border=crayon"green bold",
            ),
        )
    end
end

function pprint(SM::StoichMatrix; row_label = :symbol, col_label = :symbol, label = :identity)
    (; A, primaries, species) = SM
    pprint(A, primaries, species; row_label = row_label, col_label = col_label, label = label)
end

"""
        CanonicalStoichMatrix(species)

Construct a StoichMatrix from a list of species i.e. the matrix giving the stoichiometric
coefficients of the species in columns with respect to the atoms in rows.

# Arguments

  - `species`: list of species.

# Examples

```julia
julia> H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻ = Species.(split("H₂O H⁺ OH⁻ CO₂ HCO₃⁻ CO₃²⁻"))
6-element Vector{Species{Int64}}:
 H₂O {H₂O} [H₂O ◆ H2O]
 H⁺ {H⁺} [H⁺ ◆ H+]
 OH⁻ {OH⁻} [OH⁻ ◆ OH-]
 CO₂ {CO₂} [CO₂ ◆ CO2]
 HCO₃⁻ {HCO₃⁻} [HCO₃⁻ ◆ HCO3-]
 CO₃²⁻ {CO₃²⁻} [CO₃²⁻ ◆ CO3-2]

julia> CSM = CanonicalStoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻])
┌────┬─────┬────┬─────┬─────┬───────┬───────┐
│    │ H₂O │ H⁺ │ OH⁻ │ CO₂ │ HCO₃⁻ │ CO₃²⁻ │
├────┼─────┼────┼─────┼─────┼───────┼───────┤
│  C │     │    │     │   1 │     1 │     1 │
│  H │   2 │  1 │   1 │     │     1 │       │
│  O │   1 │    │   1 │   2 │     3 │     3 │
│ Zz │     │  1 │  -1 │     │    -1 │    -2 │
└────┴─────┴────┴─────┴─────┴───────┴───────┘
```
"""
function CanonicalStoichMatrix(species::AbstractVector{<:AbstractSpecies})
    involved_atoms_dicts = same_components(species).(species)
    involved_atoms = union_atoms(involved_atoms_dicts, item_order(species))
    T = promote_type(valtype.(involved_atoms_dicts)...)
    A = zeros(T, length(involved_atoms), length(species))
    for (j, atoms) in enumerate(involved_atoms_dicts)
        for (i, atom) in enumerate(involved_atoms)
            A[i, j] = get(atoms, atom, zero(T))
        end
    end

    return StoichMatrix(A, involved_atoms, species, similar(A, 0, 0))
end

"""
        StoichMatrix(species, candidate_primaries=species; involve_all_atoms=true)

Construct a StoichMatrix from a list of species and a list of candidate primary species
(by default the list of species itself).

# Arguments

  - `species`: list of species.
  - `candidate_primaries`: list of candidate primary species (default: `species`).
  - `involve_all_atoms`: if true the algorithm is allowed to use species
  of `candidate_primaries` containing atoms which are not in `species` (default: true).

# Examples

```julia
julia> H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻ = Species.(split("H₂O H⁺ OH⁻ CO₂ HCO₃⁻ CO₃²⁻"))
6-element Vector{Species{Int64}}:
 H₂O {H₂O} [H₂O ◆ H2O]
 H⁺ {H⁺} [H⁺ ◆ H+]
 OH⁻ {OH⁻} [OH⁻ ◆ OH-]
 CO₂ {CO₂} [CO₂ ◆ CO2]
 HCO₃⁻ {HCO₃⁻} [HCO₃⁻ ◆ HCO3-]
 CO₃²⁻ {CO₃²⁻} [CO₃²⁻ ◆ CO3-2]

julia> SM = StoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻])
┌─────┬─────┬────┬─────┬─────┬───────┬───────┐
│     │ H₂O │ H⁺ │ OH⁻ │ CO₂ │ HCO₃⁻ │ CO₃²⁻ │
├─────┼─────┼────┼─────┼─────┼───────┼───────┤
│ H₂O │   1 │    │   1 │     │     1 │     1 │
│  H⁺ │     │  1 │  -1 │     │    -1 │    -2 │
│ CO₂ │     │    │     │   1 │     1 │     1 │
└─────┴─────┴────┴─────┴─────┴───────┴───────┘

julia> SM.N # nullspace (columns are vectors of the kernel of the stoichiometric matrix)
6×3 Matrix{Int64}:
 -1  -1  -1
  1   1   2
  1   0   0
  0  -1  -1
  0   1   0
  0   0   1

julia> SM.A*SM.N
3×3 Matrix{Int64}:
 0  0  0
 0  0  0
 0  0  0
```
"""
function StoichMatrix(
    species::AbstractVector{<:AbstractSpecies},
    candidate_primaries::AbstractVector{<:AbstractSpecies}=species;
    involve_all_atoms=true,
    optimize_primaries=false
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

    all_species = union(species, candidate_primaries)
    vec_components = same_components(all_species)

    S = promote_type(typeof.(species)..., typeof.(candidate_primaries)...)

    newspecies = S[]
    append!(newspecies, species)
    num_initial_species = length(newspecies)
    initial_involved_atoms = if involve_all_atoms
        union_atoms(vec_components.(all_species), item_order(newspecies))
    else
        union_atoms(vec_components.(newspecies), item_order(newspecies))
    end
    candidate_primaries = deepcopy(candidate_primaries)

    SpType(::Vector) = Species
    SpType(::Vector{<:CemSpecies}) = CemSpecies
    Zz = SpType(newspecies)("Zz")
    charged = :Zz ∈ initial_involved_atoms
    if charged
        if Zz ∉ newspecies
            push!(newspecies, Zz)
            num_initial_species += 1
        end
        if Zz ∉ candidate_primaries
            push!(candidate_primaries, Zz)
        end
    end

    for x in candidate_primaries
        idx = findfirst(y -> x == y, newspecies)
        if isnothing(idx) &&
            all(k -> first(k) ∈ initial_involved_atoms || first(k) == :Zz, vec_components(x))
            push!(newspecies, x)
        end
    end

    CSM = CanonicalStoichMatrix(newspecies)
    M, involved_atoms = CSM.A, CSM.primaries
    redox = charged && safe_rank(M[:, begin:end .!= num_initial_species]) != safe_rank(M[1:(end - 1), begin:end .!= num_initial_species])

    if !redox && charged
        deleteat!(newspecies, num_initial_species)
        M = M[1:(end - 1), 1:end .!= num_initial_species]
        num_initial_species -= 1
    end

    cols_candidates = [findfirst(y -> y == x, newspecies) for x in candidate_primaries]
    filter!(x -> !isnothing(x), cols_candidates)
    M_subset = M[:, cols_candidates]

    r = Int(safe_rank(M_subset))
    if optimize_primaries
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
            symbol(newspecies[x]) !== "H2O@" &&
            symbol(newspecies[x]) !== "H2O" &&
            symbol(newspecies[x]) !== "H₂O" &&
            symbol(newspecies[x]) !== "H",
    )
    M_indep = M[:, independent_cols_indices]
    M_indep = promote_type(typeof.(M_indep)...).(M_indep)
    A = stoich_coef_round.(safe_pinv(M_indep) * M)

    indep_comp = newspecies[independent_cols_indices]
    dep_comp = newspecies[1:num_initial_species]
    A = A[:, 1:num_initial_species]

    if redox && Zz ∈ dep_comp && Zz ∉ indep_comp
        A = A[:, 1:(end - 1)]
        dep_comp = dep_comp[1:(end - 1)]
        num_initial_species -= 1
    end

    zero_rows = all(iszero.(A), dims=2)[:, 1]
    A = A[.!zero_rows, :]
    indep_comp = indep_comp[.!zero_rows]

    return StoichMatrix(A, indep_comp, dep_comp)
end

"""
        pull_primaries(SM::StoichMatrix)

Construct a StoichMatrix by reordering the species such that the primaries appear first
(works only if primaries are contained within species).

# Arguments

  - `SM`: StoichMatrix.

# Examples

```julia
julia> H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻ = Species.(split("H₂O H⁺ OH⁻ CO₂ HCO₃⁻ CO₃²⁻"))
6-element Vector{Species{Int64}}:
 H₂O {H₂O} [H₂O ◆ H2O]
 H⁺ {H⁺} [H⁺ ◆ H+]
 OH⁻ {OH⁻} [OH⁻ ◆ OH-]
 CO₂ {CO₂} [CO₂ ◆ CO2]
 HCO₃⁻ {HCO₃⁻} [HCO₃⁻ ◆ HCO3-]
 CO₃²⁻ {CO₃²⁻} [CO₃²⁻ ◆ CO3-2]

julia> SM = StoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻])
┌─────┬─────┬────┬─────┬─────┬───────┬───────┐
│     │ H₂O │ H⁺ │ OH⁻ │ CO₂ │ HCO₃⁻ │ CO₃²⁻ │
├─────┼─────┼────┼─────┼─────┼───────┼───────┤
│ H₂O │   1 │    │   1 │     │     1 │     1 │
│  H⁺ │     │  1 │  -1 │     │    -1 │    -2 │
│ CO₂ │     │    │     │   1 │     1 │     1 │
└─────┴─────┴────┴─────┴─────┴───────┴───────┘

julia> pull_primaries(SM)
┌─────┬─────┬────┬─────┬─────┬───────┬───────┐
│     │ H₂O │ H⁺ │ CO₂ │ OH⁻ │ HCO₃⁻ │ CO₃²⁻ │
├─────┼─────┼────┼─────┼─────┼───────┼───────┤
│ H₂O │   1 │    │     │   1 │     1 │     1 │
│  H⁺ │     │  1 │     │  -1 │    -1 │    -2 │
│ CO₂ │     │    │   1 │     │     1 │     1 │
└─────┴─────┴────┴─────┴─────┴───────┴───────┘

julia> pull_primaries(SM).N
6×3 Matrix{Int64}:
 -1  -1  -1
  1   1   2
  0  -1  -1
  1   0   0
  0   1   0
  0   0   1
```
"""
function pull_primaries(SM::StoichMatrix)
    (; A, primaries, species) = SM
    indices_in = [findfirst(x -> x == p, species) for p in primaries]
    if any(x->isnothing(x), indices_in)
        return SM
    else
        indices_out = setdiff(eachindex(species), indices_in)
        indices = [indices_in; indices_out]
        return StoichMatrix(A[:, indices], primaries, species[indices], [-A[:, indices_out]; I])
    end
end

"""
        push_primaries(SM::StoichMatrix)

Construct a StoichMatrix by reordering the species such that the primaries appear at the end
(works only if primaries are contained within species).

# Arguments

  - `SM`: StoichMatrix.

# Examples

```julia
julia> H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻ = Species.(split("H₂O H⁺ OH⁻ CO₂ HCO₃⁻ CO₃²⁻"))
6-element Vector{Species{Int64}}:
 H₂O {H₂O} [H₂O ◆ H2O]
 H⁺ {H⁺} [H⁺ ◆ H+]
 OH⁻ {OH⁻} [OH⁻ ◆ OH-]
 CO₂ {CO₂} [CO₂ ◆ CO2]
 HCO₃⁻ {HCO₃⁻} [HCO₃⁻ ◆ HCO3-]
 CO₃²⁻ {CO₃²⁻} [CO₃²⁻ ◆ CO3-2]

julia> SM = StoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻])
┌─────┬─────┬────┬─────┬─────┬───────┬───────┐
│     │ H₂O │ H⁺ │ OH⁻ │ CO₂ │ HCO₃⁻ │ CO₃²⁻ │
├─────┼─────┼────┼─────┼─────┼───────┼───────┤
│ H₂O │   1 │    │   1 │     │     1 │     1 │
│  H⁺ │     │  1 │  -1 │     │    -1 │    -2 │
│ CO₂ │     │    │     │   1 │     1 │     1 │
└─────┴─────┴────┴─────┴─────┴───────┴───────┘

julia> push_primaries(SM)
┌─────┬─────┬───────┬───────┬─────┬────┬─────┐
│     │ OH⁻ │ HCO₃⁻ │ CO₃²⁻ │ H₂O │ H⁺ │ CO₂ │
├─────┼─────┼───────┼───────┼─────┼────┼─────┤
│ H₂O │   1 │     1 │     1 │   1 │    │     │
│  H⁺ │  -1 │    -1 │    -2 │     │  1 │     │
│ CO₂ │     │     1 │     1 │     │    │   1 │
└─────┴─────┴───────┴───────┴─────┴────┴─────┘

julia> push_primaries(SM).N
6×3 Matrix{Int64}:
  1   0   0
  0   1   0
  0   0   1
 -1  -1  -1
  1   1   2
  0  -1  -1
```
"""
function push_primaries(SM::StoichMatrix)
    (; A, primaries, species) = SM
    indices_in = [findfirst(x -> x == p, species) for p in primaries]
    if any(x->isnothing(x), indices_in)
        return SM
    else
        indices_out = setdiff(eachindex(species), indices_in)
        indices = [indices_out; indices_in]
        return StoichMatrix(A[:, indices], primaries, species[indices], [I; -A[:, indices_out]])
    end
end

"""
        mass_matrix(SM::StoichMatrix)

Construct a StoichMatrix with mass correspondence instead of molar stoichiometry.

# Arguments

  - `SM`: StoichMatrix.

# Examples

```julia
julia> H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻ = Species.(split("H₂O H⁺ OH⁻ CO₂ HCO₃⁻ CO₃²⁻"))
6-element Vector{Species{Int64}}:
 H₂O {H₂O} [H₂O ◆ H2O]
 H⁺ {H⁺} [H⁺ ◆ H+]
 OH⁻ {OH⁻} [OH⁻ ◆ OH-]
 CO₂ {CO₂} [CO₂ ◆ CO2]
 HCO₃⁻ {HCO₃⁻} [HCO₃⁻ ◆ HCO3-]
 CO₃²⁻ {CO₃²⁻} [CO₃²⁻ ◆ CO3-2]

julia> SM = StoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻])
┌─────┬─────┬────┬─────┬─────┬───────┬───────┐
│     │ H₂O │ H⁺ │ OH⁻ │ CO₂ │ HCO₃⁻ │ CO₃²⁻ │
├─────┼─────┼────┼─────┼─────┼───────┼───────┤
│ H₂O │   1 │    │   1 │     │     1 │     1 │
│  H⁺ │     │  1 │  -1 │     │    -1 │    -2 │
│ CO₂ │     │    │     │   1 │     1 │     1 │
└─────┴─────┴────┴─────┴─────┴───────┴───────┘

julia> mass_matrix(SM)
┌─────┬─────┬─────┬────────────┬─────┬────────────┬────────────┐
│     │ H₂O │  H⁺ │        OH⁻ │ CO₂ │      HCO₃⁻ │      CO₃²⁻ │
├─────┼─────┼─────┼────────────┼─────┼────────────┼────────────┤
│ H₂O │ 1.0 │     │    1.05927 │     │    0.29525 │    0.30021 │
│  H⁺ │     │ 1.0 │ -0.0592697 │     │ -0.0165203 │ -0.0335955 │
│ CO₂ │     │     │            │ 1.0 │    0.72127 │   0.733386 │
└─────┴─────┴─────┴────────────┴─────┴────────────┴────────────┘

julia> CSM = CanonicalStoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻])
┌────┬─────┬────┬─────┬─────┬───────┬───────┐
│    │ H₂O │ H⁺ │ OH⁻ │ CO₂ │ HCO₃⁻ │ CO₃²⁻ │
├────┼─────┼────┼─────┼─────┼───────┼───────┤
│  C │     │    │     │   1 │     1 │     1 │
│  H │   2 │  1 │   1 │     │     1 │       │
│  O │   1 │    │   1 │   2 │     3 │     3 │
│ Zz │     │  1 │  -1 │     │    -1 │    -2 │
└────┴─────┴────┴─────┴─────┴───────┴───────┘

julia> mass_matrix(CSM)
┌────┬──────────┬─────┬───────────┬──────────┬───────────┬──────────┐
│    │      H₂O │  H⁺ │       OH⁻ │      CO₂ │     HCO₃⁻ │    CO₃²⁻ │
├────┼──────────┼─────┼───────────┼──────────┼───────────┼──────────┤
│  C │          │     │           │ 0.272921 │   0.19685 │ 0.200157 │
│  H │ 0.111907 │ 1.0 │ 0.0592697 │          │ 0.0165203 │          │
│  O │ 0.888093 │     │   0.94073 │ 0.727079 │   0.78663 │ 0.799843 │
│ Zz │          │     │           │          │           │          │
└────┴──────────┴─────┴───────────┴──────────┴───────────┴──────────┘
```
"""
function mass_matrix(SM::StoichMatrix{T,S}) where {T,S<:AbstractSpecies}
    Mspecies = ustrip.(getproperty.(SM.species, :M))
    Mprimaries = ustrip.(getproperty.(SM.primaries, :M))
    return StoichMatrix(Mprimaries .* SM.A .* inv.(Mspecies)', SM.primaries, SM.species)
end

function mass_matrix(SM::StoichMatrix{T,Symbol}) where {T}
    Mspecies = ustrip.(getproperty.(SM.species, :M))
    S = promote_type(root_type.(typeof.(SM.species))...)
    Mprimaries = ustrip.(getproperty.(S.(string.(SM.primaries)), :M))
    return StoichMatrix(Mprimaries .* SM.A .* inv.(Mspecies)', SM.primaries, SM.species)
end

"""
        reactions(SM::StoichMatrix)

Construct the list of non trivial reactions from a StoichMatrix.

# Arguments

  - `SM`: StoichMatrix.

# Examples

```julia
julia> H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻ = Species.(split("H₂O H⁺ OH⁻ CO₂ HCO₃⁻ CO₃²⁻"))
6-element Vector{Species{Int64}}:
 H₂O {H₂O} [H₂O ◆ H2O]
 H⁺ {H⁺} [H⁺ ◆ H+]
 OH⁻ {OH⁻} [OH⁻ ◆ OH-]
 CO₂ {CO₂} [CO₂ ◆ CO2]
 HCO₃⁻ {HCO₃⁻} [HCO₃⁻ ◆ HCO3-]
 CO₃²⁻ {CO₃²⁻} [CO₃²⁻ ◆ CO3-2]

julia> SM = StoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻])
┌─────┬─────┬────┬─────┬─────┬───────┬───────┐
│     │ H₂O │ H⁺ │ OH⁻ │ CO₂ │ HCO₃⁻ │ CO₃²⁻ │
├─────┼─────┼────┼─────┼─────┼───────┼───────┤
│ H₂O │   1 │    │   1 │     │     1 │     1 │
│  H⁺ │     │  1 │  -1 │     │    -1 │    -2 │
│ CO₂ │     │    │     │   1 │     1 │     1 │
└─────┴─────┴────┴─────┴─────┴───────┴───────┘

julia> reactions(SM)
3-element Vector{Reaction{Species{Int64}, Int64, Species{Int64}, Int64, Int64}}:
 H₂O = OH⁻ + H⁺
 H₂O + CO₂ = HCO₃⁻ + H⁺
 H₂O + CO₂ = CO₃²⁻ + 2H⁺
```
"""
function reactions(SM::StoichMatrix)
    if !isempty(SM.N)
        pSM = push_primaries(SM)
        return [Reaction(OrderedDict(zip(pSM.species, V)); symbol=symbol(pSM.species[j])) for (j, V) in enumerate(eachcol(pSM.N))]
    else
        lr = unique!([Reaction(merge(+, OrderedDict(SM.species[j]=>1), OrderedDict(zip(SM.primaries, -V))); symbol=symbol(SM.species[j])) for (j, V) in enumerate(eachcol(SM.A))])
        return filter(x->!isempty(x), lr)
    end
end

"""
        pprint(reactions::AbstractVector{<:Reaction})

Pretty print a list of reactions.

# Arguments

  - `SM`: StoichMatrix.

# Examples

```julia
julia> H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻ = Species.(split("H₂O H⁺ OH⁻ CO₂ HCO₃⁻ CO₃²⁻"))
6-element Vector{Species{Int64}}:
 H₂O {H₂O} [H₂O ◆ H2O]
 H⁺ {H⁺} [H⁺ ◆ H+]
 OH⁻ {OH⁻} [OH⁻ ◆ OH-]
 CO₂ {CO₂} [CO₂ ◆ CO2]
 HCO₃⁻ {HCO₃⁻} [HCO₃⁻ ◆ HCO3-]
 CO₃²⁻ {CO₃²⁻} [CO₃²⁻ ◆ CO3-2]

julia> SM = StoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻])
┌─────┬─────┬────┬─────┬─────┬───────┬───────┐
│     │ H₂O │ H⁺ │ OH⁻ │ CO₂ │ HCO₃⁻ │ CO₃²⁻ │
├─────┼─────┼────┼─────┼─────┼───────┼───────┤
│ H₂O │   1 │    │   1 │     │     1 │     1 │
│  H⁺ │     │  1 │  -1 │     │    -1 │    -2 │
│ CO₂ │     │    │     │   1 │     1 │     1 │
└─────┴─────┴────┴─────┴─────┴───────┴───────┘

julia> pprint(reactions(SM))
  OH⁻ │ H₂O = OH⁻ + H⁺
HCO₃⁻ │ H₂O + CO₂ = HCO₃⁻ + H⁺
CO₃²⁻ │ H₂O + CO₂ = CO₃²⁻ + 2H⁺
```
"""
function pprint(reactions::AbstractVector{<:Reaction})
    pad = maximum(length.(symbol.(reactions)))
    for r in reactions
        println(lpad("$(symbol(r))", pad), " │ ", colored(r))
    end
end

const oxides_as_species = [Species(d; symbol=string(k)) for (k, d) in CEMENT_TO_MENDELEEV]

const Aoxides, atoms_in_oxides = getfield.(Ref(CanonicalStoichMatrix(oxides_as_species)), [:A, :primaries])

const order_atom_in_oxides = Dict(atom => i for (i, atom) in enumerate(atoms_in_oxides))
