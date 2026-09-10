# The page tree, grouped by what a chapter is *about* rather than by the order the
# pages were written in.
#
# Four chapters carry the documentation, and the difference between them is the
# question each answers.
#
#   **Theory**    — why is this the right calculation? Readable without running
#                   anything; every equation is the one the code evaluates.
#   **Manual**    — how is this written? One page per kind of object, syntax
#                   first, on the smallest example that shows it.
#   **Tutorials** — how do I drive a calculation from end to end? Narrative, and
#                   each one goes somewhere.
#   **Applications** — what does a real case look like, and what do the choices
#                   cost in numbers? Everything executed lives here.
#   **API**       — the docstrings, generated.
#
# A page that answers two of those questions is worth splitting; a page in the
# wrong chapter is worth moving. `Manual` exists because nine of the fourteen
# pages once filed under `Tutorials` were object syntax rather than a
# calculation, which is what made the chapter hard to navigate.
#
# **Cementitious media is a subsection of four chapters, not a chapter of its
# own.** The material is the subject of the package, so it appears wherever its
# question is being asked: the theory of its water budget, the syntax of its
# species, the trajectory of its hydration, the worked cases. Grouping it into
# one chapter would have separated each of those from the general treatment it
# specializes.
#
# Three paths are load-bearing and must not move: `tutorials/self_desiccation.md`,
# `tutorials/equilibrium.md` and `examples/hydration_calibration.md` are linked by
# URL from released CHANGELOG sections, which are not retro-edited.

pages = [
    "Home" => "index.md",
    "Getting Started" => "quickstart.md",
    "Theory" => [
        "theory/index.md",
        # The definitions and identities the rest is written in. Read first: the
        # remaining pages use its notation, which is the code's.
        "Foundations" => [
            "theory/thermodynamics.md",
            "theory/equilibrium.md",
        ],
        # The two places a mixture stops being ideal, and the only two where a
        # standard state has to be argued about rather than looked up.
        "Non-ideal mixtures" => [
            "theory/activity_models.md",
            "theory/solid_solutions.md",
        ],
        # The material this package exists for. What a Gibbs minimization can
        # predict about a drying paste, and what is not a thermodynamic
        # quantity at all.
        "Cementitious media" => [
            "theory/cement_water_budget.md",
        ],
    ],
    "Manual" => [
        # What a formula, a species and a reaction *are* here. Everything else
        # consumes these.
        "Chemical description" => [
            "manual/formula_manipulation.md",
            "manual/species.md",
            "manual/reactions.md",
            "manual/stoich_matrices.md",
        ],
        "Databases and thermodynamic data" => [
            "manual/databases.md",
            "manual/thermodynamic_data.md",
        ],
        "Systems and states" => [
            "manual/chemical_system_state.md",
        ],
        # The cement-specific syntax: phase names, the Bogue notation, the
        # shorthand the literature uses.
        "Cementitious media" => [
            "manual/cement_species.md",
        ],
        "Appendices" => [
            "manual/advanced.md",
        ],
    ],
    "Tutorials" => [
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
        "Validation against other codes" => [
            "tutorials/reaktoro_comparison.md",
        ],
    ],
    "Applications" => [
        # From an oxide analysis to a species list — the entry point for someone
        # holding a cement datasheet rather than a database.
        "From a cement analysis to a chemical system" => [
            "examples/bogue_calculation.md",
            "examples/example_stoich_matrix.md",
            "examples/from_scratch.md",
        ],
        # The executed counterpart of the Theory chapter: what the modeling
        # choices cost, in numbers, on a composition small enough to check.
        "Non-ideal mixtures, measured" => [
            "examples/activity_models_compared.md",
            "examples/solid_solution_models.md",
        ],
        # Small, checkable aqueous cases with an analytical answer to compare to.
        "Aqueous equilibria" => [
            "examples/titration_acetic_acid.md",
            "examples/titration_malonic_acid.md",
            "examples/co2_carbonate_system.md",
        ],
        "Cementitious media at equilibrium" => [
            "examples/simplified_clinker_dissolution.md",
            "examples/cement_wc_ratio.md",
            "examples/cement_carbonation.md",
        ],
        "Cementitious media in time" => [
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
