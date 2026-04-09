__precompile__(true)

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

# Examples

```julia
julia> using ChemistryLab

julia> f = Formula("H2O")
Formula{Int64}
    formula: H2O ◆ H₂O ◆ H₂O
composition: H => 2, O => 1
     charge: 0

julia> f[:H]
2

julia> H2O = Species(f; name="Vapour", symbol="H₂O", aggregate_state=AS_GAS, class=SC_GASFLUID)
Species{Int64}
           name: Vapour
         symbol: H₂O
        formula: H2O ◆ H₂O ◆ H₂O
          atoms: H => 2, O => 1
         charge: 0
aggregate_state: AS_GAS
          class: SC_GASFLUID
     properties: M = 0.0180149999937744 kg mol⁻¹

julia> CO2 = Species(
           Dict(:C=>1, :O=>2);
           name="Carbon dioxide",
           symbol="CO₂",
           aggregate_state=AS_GAS,
           class=SC_GASFLUID,
       ) # definition from Dict
Species{Int64}
           name: Dioxygen
         symbol: CO₂
        formula: CO2 ◆ CO₂ ◆ CO₂
          atoms: O => 2, C => 1
         charge: 0
aggregate_state: AS_GAS
          class: SC_GASFLUID
     properties: M = 0.04400899998479143 kg mol⁻¹

julia> O2 = Species("O2"; name="Dioxygen", symbol="O₂", aggregate_state=AS_GAS, class=SC_GASFLUID) # definition from String
Species{Int64}
           name: Carbon dioxide
         symbol: CO₂
        formula: CO2 ◆ CO₂ ◆ CO₂
          atoms: C => 1, O => 2
         charge: 0
aggregate_state: AS_GAS
          class: SC_GASFLUID
     properties: M = 0.04400899998479143 kg mol⁻¹

julia> C3H8 = Species(
           "C₃H₈"; name="Propane", symbol="C₃H₈", aggregate_state=AS_GAS, class=SC_GASFLUID
       ) # definition from Unicode String
Species{Int64}
           name: Propane
         symbol: C₃H₈
        formula: C₃H₈ ◆ C3H8 ◆ C₃H₈
          atoms: C => 3, H => 8
         charge: 0
aggregate_state: AS_GAS
          class: SC_GASFLUID
     properties: M = 0.04409699998476102 kg mol⁻¹

julia> r = Reaction([C3H8, O2, CO2, H2O])
C₃H₈ + 5O₂ = 4H₂O + 3CO₂
 reactants: C₃H₈ => 1, O₂ => 5
  products: H₂O => 4, CO₂ => 3
```
"""
module ChemistryLab

    import Base: ==, +, -, *, /, //, ^

    using Crayons
    using DataFrames
    using DynamicQuantities
    using ForwardDiff
    using JSON
    using LinearAlgebra
    using Symbolics
    using OrderedCollections
    using Tables
    using PeriodicTable
    using PrettyTables
    using ProgressMeter
    using RuntimeGeneratedFunctions
    using SciMLBase
    using TOML
    using Unicode

    RuntimeGeneratedFunctions.init(@__MODULE__)

    include("utils/misc.jl")
    include("utils/subsuperscripts.jl")

    include("thermodynamics/thermo_factories.jl")
    include("thermodynamics/water_properties.jl")
    include("thermodynamics/thermo_models.jl")
    include("thermodynamics/precompile.jl")

    include("chemical_structs/element_order.jl")
    include("chemical_structs/parsing_tools.jl")
    include("chemical_structs/formulas.jl")
    include("chemical_structs/species.jl")
    include("chemical_structs/solid_solutions.jl")
    include("chemical_structs/reactions.jl")
    include("chemical_structs/speciation.jl")
    include("chemical_structs/stoich_matrices.jl")
    include("chemical_structs/chemical_systems.jl")
    include("chemical_structs/chemical_states.jl")

    include("databases/phreeqc_dat.jl")
    include("databases/thermofun_json.jl")
    include("databases/merge_dat_json.jl")

    include("equilibrium/activities.jl")
    include("equilibrium/equilibrium_problems.jl")
    include("equilibrium/equilibrium_solver.jl")

    export SymbolicFunc,
        ThermoFactory,
        NumericFunc,
        AbstractFunc,
        infer_unit,
        derivative

    export WaterThermoProps,
        WaterElectroProps,
        HKFGState,
        SpeciesElectroPropsHKF,
        water_thermo_props,
        water_electro_props_jn,
        hkf_g_function,
        species_electro_props_hkf

    export THERMO_MODELS,
        THERMO_FACTORIES,
        add_thermo_model,
        build_thermo_functions,
        check_dimensions

    export ATOMIC_ORDER,
        CEMENT_TO_MENDELEEV,
        OXIDE_ORDER,
        CEMDATA_PRIMARIES

    export stoich_coef_round,
        phreeqc_to_unicode,
        unicode_to_phreeqc,
        colored_formula,
        parse_formula,
        extract_charge,
        to_mendeleev,
        parse_equation,
        colored_equation,
        format_equation

    export AtomGroup,
        Formula,
        expr,
        phreeqc,
        unicode,
        colored,
        composition,
        charge,
        check_mendeleev,
        calculate_molar_mass,
        stoichtype,
        pprint

    export AggregateState,
        AS_UNDEF,
        AS_AQUEOUS,
        AS_CRYSTAL,
        AS_GAS

    export Class,
        SC_UNDEF,
        SC_AQSOLVENT,
        SC_AQSOLUTE,
        SC_COMPONENT,
        SC_GASFLUID,
        SC_SSENDMEMBER

    export AbstractSpecies,
        Species,
        CemSpecies,
        name,
        symbol,
        formula,
        cemformula,
        atoms,
        atoms_charge,
        oxides,
        oxides_charge,
        components,
        aggregate_state,
        class,
        properties,
        apply

    export AbstractReaction,
        Reaction,
        CemReaction,
        reactants,
        products,
        charge,
        simplify_reaction
    @eval export $(Symbol.(EQUAL_OPS)...)

    export union_atoms, speciation

    export StoichMatrix,
        CanonicalStoichMatrix,
        pull_primaries,
        push_primaries,
        mass_matrix,
        reactions

    export AbstractSolidSolutionModel,
        IdealSolidSolutionModel,
        RedlichKisterModel,
        AbstractSolidSolutionPhase,
        SolidSolutionPhase,
        end_members,
        model,
        with_class

    export ChemicalSystem,
        aqueous,
        crystal,
        gas,
        solutes,
        solvent,
        components,
        gasfluid,
        get_reaction,
        solid_solutions

    export ChemicalState,
        PhaseQuantities,
        temperature,
        pressure,
        set_temperature!,
        set_pressure!,
        moles,
        set_quantity!,
        rescale!,
        mass,
        volume,
        pH,
        pOH,
        porosity,
        saturation

    export extract_primary_species

    export read_thermofun_database,
        build_species,
        build_reactions,
        build_solid_solutions,
        get_compatible_species,
        HKF_SI_CONVERSIONS

    export merge_json

    export AbstractActivityModel,
        DiluteSolutionModel,
        HKFActivityModel,
        DaviesActivityModel,
        activity_model,
        build_potentials,
        REJ_HKF,
        REJ_CHARGE_DEFAULT,
        hkf_debye_huckel_params

    export EquilibriumProblem

    export EquilibriumSolver,
        equilibrate

    function __init__()
        for (k, v) in THERMO_MODELS
            THERMO_FACTORIES[k] = build_thermo_factories(v)
        end
        return
    end

end
