using TOML

@testsection "Databases" begin
    # Test parse_reaction_stoich_cemdata
    @testset "parse_reaction_stoich_cemdata" begin
        # Test basic reaction parsing
        reaction = "CaCO3 = Ca+2 + CO3-2"
        reactants, equation, comment = ChemistryLab.parse_reaction_stoich_cemdata(reaction)
        @test length(reactants) == 3
        @test any(r -> r["symbol"] == "CaCO3" && r["coefficient"] == -1.0, reactants)
        @test any(r -> r["symbol"] == "Ca+2" && r["coefficient"] == 1.0, reactants)
        @test any(r -> r["symbol"] == "CO3-2" && r["coefficient"] == 1.0, reactants)

        # Test reaction with comment
        reaction_with_comment = "H2O = H+ + OH- # water dissociation"
        reactants, equation, comment = ChemistryLab.parse_reaction_stoich_cemdata(reaction_with_comment)
        @test comment == "water dissociation"
        @test length(reactants) == 3

        # Test reaction with coefficients
        reaction_with_coef = "2H2O = 2H+ + 2OH-"
        reactants, equation, comment = ChemistryLab.parse_reaction_stoich_cemdata(reaction_with_coef)
        @test length(reactants) == 3
        @test any(r -> r["symbol"] == "H2O" && r["coefficient"] == -2.0, reactants)
    end

    # Test parse_float_array
    @testset "parse_float_array" begin
        line = "-analytical_expression 1.23 -4.56 7.89 # some comment"
        values = ChemistryLab.parse_float_array(line)
        @test length(values) == 3
        @test values ≈ [1.23, -4.56, 7.89]

        # Test empty line
        @test isempty(ChemistryLab.parse_float_array(""))

        # Test line with only comments
        @test isempty(ChemistryLab.parse_float_array("# only comment"))
    end

    # Test parse_phases
    @testset "parse_phases" begin
        dat_content = """
        PHASES
        Calcite
        CaCO3 = Ca+2 + CO3-2
        -log_K -8.48
        -analytical_expression 1.23 -4.56 7.89

        Portlandite
        Ca(OH)2 = Ca+2 + 2OH-
        -log_K -5.2
        """

        phases = ChemistryLab.parse_phases(dat_content)
        @test haskey(phases, "Calcite")
        @test haskey(phases, "Portlandite")
        @test phases["Calcite"]["logKr"]["values"][1] ≈ -8.48
        @test length(phases["Calcite"]["analytical_expression"]) == 3
    end
end

@testsection "build_solid_solutions" begin
    # ── Helpers: build a minimal species dict without loading real databases ──
    _make(sym, cls = SC_COMPONENT) =
        Species(sym; aggregate_state = AS_CRYSTAL, class = cls)

    # Build a fake dict that mimics what build_species returns
    dict = Dict(
        "TobD" => _make("TobD"),
        "TobH" => _make("TobH"),
        "JenH" => _make("JenH"),
        "JenD" => _make("JenD"),
        "Ms"   => _make("Ms"),
        "Mc"   => _make("Mc"),
    )

    # ── Write a temp TOML ────────────────────────────────────────────────────
    toml_content = """
    [[solid_solution]]
    name        = "CSHQ"
    end_members = ["TobD", "TobH", "JenH", "JenD"]
    model       = "ideal"

    [[solid_solution]]
    name        = "AFm"
    end_members = ["Ms", "Mc"]
    model       = "redlich_kister"
    a0          = 3000.0
    a1          = 500.0
    a2          = 0.0

    [[solid_solution]]
    name        = "Missing_phase"
    end_members = ["NonExistent1", "NonExistent2"]
    model       = "ideal"
    """
    tmp = tempname() * ".toml"
    write(tmp, toml_content)

    @testset "load two phases, skip missing" begin
        phases = build_solid_solutions(tmp, dict; skip_missing = true)
        @test length(phases) == 2       # Missing_phase skipped

        cshq = first(filter(ss -> name(ss) == "CSHQ", phases))
        afm  = first(filter(ss -> name(ss) == "AFm",  phases))

        @test length(end_members(cshq)) == 4
        @test model(cshq) isa IdealSolidSolutionModel

        @test length(end_members(afm)) == 2
        @test model(afm) isa RedlichKisterModel
        @test model(afm).a0 ≈ 3000.0
        @test model(afm).a1 ≈ 500.0
    end

    @testset "end-members requalified to SC_SSENDMEMBER" begin
        phases = build_solid_solutions(tmp, dict; skip_missing = true)
        for ss in phases
            for em in end_members(ss)
                @test class(em) == SC_SSENDMEMBER
            end
        end
    end

    @testset "original dict species class untouched" begin
        build_solid_solutions(tmp, dict; skip_missing = true)
        # with_class returns a copy; dict values must still be SC_COMPONENT
        @test all(class(s) == SC_COMPONENT for s in values(dict))
    end

    @testset "skip_missing = false raises on missing end-members" begin
        @test_throws ErrorException build_solid_solutions(tmp, dict; skip_missing = false)
    end

    @testset "data/solid_solutions.toml parses without error" begin
        toml_path = joinpath(@__DIR__, "..", "data", "solid_solutions.toml")
        # Just check the file is valid TOML and has solid_solution entries
        data = TOML.parsefile(toml_path)
        @test haskey(data, "solid_solution")
        @test length(data["solid_solution"]) >= 2
        @test any(e -> e["name"] == "CSHQ", data["solid_solution"])
        @test any(e -> e["name"] == "AFm",  data["solid_solution"])
    end

    rm(tmp; force = true)
end
