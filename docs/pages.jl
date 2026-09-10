# The page tree, grouped by what a chapter is *about* rather than by the order
# the pages were written in.
#
# Three chapters carry the documentation, and the difference between them is the
# question each answers. **Theory** answers "why is this the right calculation",
# and its pages are readable without running anything. **Tutorials** answer "how
# do I drive this feature", one feature at a time, on the smallest system that
# shows it. **Examples** answer "what does a real case look like end to end", and
# they are long on purpose. A page that started in the wrong chapter is worth
# moving; a page that answers two of those questions is worth splitting.
#
# Within a chapter the subsections follow the dependency chain — formulas and
# species before the systems built from them, equilibrium before the kinetics
# that calls it, and the cement material last, because it uses all of it.

pages = [
    "Home" => "index.md",
    "Getting Started" => "quickstart.md",
    "Theory" => [
        "theory/index.md",
        # The definitions and identities everything else is written in. Read
        # first: the rest of the chapter uses its notation, which is the code's.
        "Foundations" => [
            "theory/thermodynamics.md",
        ],
        # The two places a mixture stops being ideal, and the only two places a
        # standard state has to be argued about rather than looked up.
        "Non-ideal mixtures" => [
            "theory/activity_models.md",
            "theory/solid_solutions.md",
        ],
    ],
    "Tutorials" => [
        # What a species and a reaction *are* in this package, and where they
        # come from. Everything else consumes these.
        "Chemical description" => [
            "tutorials/formula_manipulation.md",
            "tutorials/species.md",
            "tutorials/cement_species.md",
            "tutorials/databases.md",
            "tutorials/stoich_matrices.md",
            "tutorials/reactions.md",
        ],
        # Standard-state data, then the containers the solvers act on.
        "Thermodynamic data, systems and states" => [
            "tutorials/thermodynamics.md",
            "tutorials/chemical_system_state.md",
        ],
        "Equilibrium" => [
            "tutorials/equilibrium.md",
        ],
        # The kinetics calls the equilibrium solver, so it reads after it.
        "Kinetics and coupling" => [
            "tutorials/kinetics.md",
            "tutorials/coupling.md",
        ],
        # Uses the whole chain, which is why it comes last.
        "Cementitious media" => [
            "tutorials/self_desiccation.md",
        ],
        "Comparisons and advanced use" => [
            "tutorials/reaktoro_comparison.md",
            "tutorials/advanced.md",
        ],
    ],
    "Examples" => [
        # From an oxide analysis to a species list — the entry point for someone
        # holding a cement datasheet rather than a database.
        "From a cement analysis to a chemical system" => [
            "examples/bogue_calculation.md",
            "examples/example_stoich_matrix.md",
            "examples/from_scratch.md",
        ],
        # Small, checkable aqueous cases with an analytical answer to compare to.
        "Aqueous equilibria" => [
            "examples/titration_acetic_acid.md",
            "examples/titration_malonic_acid.md",
            "examples/co2_carbonate_system.md",
        ],
        "Hydration at equilibrium" => [
            "examples/simplified_clinker_dissolution.md",
            "examples/cement_wc_ratio.md",
            "examples/cement_carbonation.md",
        ],
        "Hydration in time" => [
            "examples/cement_clinker_kinetics.md",
            "examples/coupled_hydration.md",
            "examples/ionic_hydration.md",
            "examples/hydration_calibration.md",
        ],
    ],
    "API" => Any[
        "Chemical description" => [
            "Formulas" => "api/formulas.md",
            "Species" => "api/species.md",
            "Parsing tools" => "api/parsing_tools.md",
            "Element order" => "api/element_order.md",
            "Stoichiometric Matrix" => "api/stoich_matrices.md",
            "Reactions" => "api/reactions.md",
            "Databases" => "api/databases.md",
        ],
        "Thermodynamics" => [
            "Thermodynamical functions" => "api/thermo_functions.md",
            "Thermodynamical models" => "api/thermo_models.md",
            "Water properties" => "api/water_properties.md",
        ],
        "Systems, states and equilibrium" => [
            "Chemical systems and states" => "api/chemical_systems.md",
            "Equilibrium" => "api/equilibrium.md",
        ],
        "Kinetics" => [
            "Kinetics" => "api/kinetics.md",
        ],
        "Utilities" => [
            "Utilities" => "api/utils.md",
        ],
    ],
    "References" => "references.md",
]
