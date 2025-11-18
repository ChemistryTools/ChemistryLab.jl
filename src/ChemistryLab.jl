"""
    ChemistryLab

Top-level module for parsing, representing and manipulating chemical
formulas, species, stoichiometric matrices and ThermoFun / PHREEQC-like data.

# Overview

  - Canonical containers: `AtomGroup`, `Formula`.
  - Species representations: `Species`, `CemSpecies`.
  - Parsers and IO helpers for ThermoFun / PHREEQC data (substances, reactions).
  - Stoichiometric matrix construction and reaction helpers.
  - Utilities for units, colored and Unicode terminal output.

# Note

This file mainly aggregates and re-exports functionality implemented
in submodules (parsing_utils.jl, formulas.jl, species.jl, databases.jl, ...).
Consult those files for detailed docs of individual types and functions.

# Examples

```jldoctest
julia> using ChemistryLab

julia> f = Formula("H2O");

julia> f[:H]
2
```
"""
module ChemistryLab

import Base: ==, +, -, *, /, //, ^

using Crayons
using CSV
using DataFrames
using DynamicQuantities
using JSON
using JSON3
using LinearAlgebra
using ModelingToolkit
using OrderedCollections
using PeriodicTable
using PrettyTables
using ProgressMeter
using SymbolicNumericIntegration
using Unicode

include("parsing_utils.jl")
include("thermodynamical_functions.jl")
include("formulas.jl")
include("species.jl")
include("databases.jl")
include("stoich_matrix.jl")
include("reactions.jl")

export ATOMIC_ORDER,
    OXIDE_ORDER,
    stoich_coef_round,
    phreeqc_to_unicode,
    colored_formula,
    unicode_to_phreeqc,
    parse_formula,
    extract_charge,
    calculate_molar_mass,
    to_mendeleev,
    parse_equation,
    format_equation,
    colored_equation
export Callable, ThermoFunction, ∂, ∫
export AtomGroup,
    Formula, expr, phreeqc, unicode, colored, composition, charge, apply, check_mendeleev
export AggregateState, AS_UNDEF, AS_AQUEOUS, AS_CRYSTAL, AS_GAS
export Class, SC_UNDEF, SC_AQSOLVENT, SC_AQSOLUTE, SC_COMPONENT, SC_GASFLUID
export AbstractSpecies, Species, CemSpecies
export name,
    symbol,
    formula,
    mendeleev_filter,
    cemformula,
    atoms,
    atoms_charge,
    oxides,
    oxides_charge,
    components,
    aggregate_state,
    class,
    properties
export merge_json,
    read_thermofun_substances,
    read_thermofun_reactions,
    read_thermofun_elements,
    read_thermofun,
    extract_primary_species
export union_atoms,
    print_stoich_matrix,
    canonical_stoich_matrix,
    stoich_matrix,
    stoich_matrix_to_equations,
    stoich_matrix_to_reactions
export Reaction, CemReaction, reactants, products, simplify_reaction
@eval export $(Symbol.(EQUAL_OPS)...)

end
