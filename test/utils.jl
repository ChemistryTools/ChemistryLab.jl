@testsection "Utils" begin
    @testsection "safe_ustrip" begin
        # Quantity → strip with target unit
        @test ChemistryLab.safe_ustrip(1u"m", 5u"m") ≈ 5.0
        @test ChemistryLab.safe_ustrip(1u"cm", 1u"m") ≈ 100.0

        # Plain number → returned unchanged
        @test ChemistryLab.safe_ustrip(1u"m", 3.0) ≈ 3.0
        @test ChemistryLab.safe_ustrip(1u"K", 42) == 42
    end

    @testsection "safe_uconvert" begin
        result = ChemistryLab.safe_uconvert(us"m", 5u"cm")
        @test isapprox(ustrip(result), 0.05; rtol = 1.0e-10)

        # Plain number → returned unchanged
        @test ChemistryLab.safe_uconvert(us"m", 3.0) == 3.0
    end

    @testsection "safe_uparse" begin
        # String input → parsed unit quantity
        @test ChemistryLab.safe_uparse("m/s") isa AbstractQuantity

        # Quantity input → returned as-is
        q = 1u"m"
        @test ChemistryLab.safe_uparse(q) === q
    end

    @testsection "force_uconvert" begin
        result = ChemistryLab.force_uconvert(1u"m", 100u"cm")
        @test ustrip(result) ≈ 1.0
        @test dimension(result) == dimension(u"m")

        # Plain number scaled to unit
        result2 = ChemistryLab.force_uconvert(1u"K", 300.0)
        @test ustrip(result2) ≈ 300.0
    end

    @testsection "root_type" begin
        # Vector{Int} is an alias for Array{Int,1}, so root_type returns Array
        @test ChemistryLab.root_type(Vector{Int}) == Array
        @test ChemistryLab.root_type(Dict{Symbol, Int}) == Dict
    end

    @testsection "print_title does not error" begin
        @test_nowarn ChemistryLab.print_title("Test"; style = :none)
        @test_nowarn ChemistryLab.print_title("Test"; style = :underline)
        @test_nowarn ChemistryLab.print_title("Test"; style = :box)
    end
end
