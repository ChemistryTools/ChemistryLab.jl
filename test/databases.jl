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
