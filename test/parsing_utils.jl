using ChemistryLab
using Test

@testset "Parsing Utils" begin
    @testset "Character conversion" begin
        # Superscript conversions
        @test ChemistryLab.super_to_normal("H₂O²⁺") == "H₂O2+"
        @test ChemistryLab.normal_to_super("2+") == "²⁺"

        # Subscript conversions
        @test ChemistryLab.sub_to_normal("H₂O") == "H2O"
        @test ChemistryLab.normal_to_sub("H2O") == "H₂O"

        # all_normal_to_sub
        @test ChemistryLab.all_normal_to_sub("H2O") == "H₂O"
        @test ChemistryLab.all_normal_to_sub("SO4") == "SO₄"

        # issuperscript / issubscript
        @test ChemistryLab.issuperscript('²')
        @test ChemistryLab.issuperscript('⁺')
        @test !ChemistryLab.issuperscript('a')
        @test ChemistryLab.issubscript('₂')
        @test !ChemistryLab.issubscript('2')

        # subscriptnumber / superscriptnumber round-trips
        @test ChemistryLab.subscriptnumber(42) == "₄₂"
        @test ChemistryLab.subscriptnumber(-3) == "₋₃"
        @test ChemistryLab.superscriptnumber(42) == "⁴²"
        @test ChemistryLab.superscriptnumber(-2) == "⁻²"

        # from_subscriptnumber / from_superscriptnumber
        @test ChemistryLab.from_subscriptnumber("₄₂") == 42
        @test ChemistryLab.from_subscriptnumber("₋₃") == -3
        @test ChemistryLab.from_superscriptnumber("⁴²") == 42
        @test ChemistryLab.from_superscriptnumber("⁻²") == -2

        # Round-trip for subscript numbers
        for i in [0, 1, 7, 12, 99, -5]
            @test ChemistryLab.from_subscriptnumber(ChemistryLab.subscriptnumber(i)) == i
        end

        # Round-trip for superscript numbers
        for i in [0, 1, 3, 10, 42, -1, -10]
            @test ChemistryLab.from_superscriptnumber(ChemistryLab.superscriptnumber(i)) == i
        end
    end

    @testset "Formula parsing" begin
        @test parse_formula("H2O") == Dict(:H => 2, :O => 1)
        @test parse_formula("Ca(OH)2") == Dict(:Ca => 1, :H => 2, :O => 2)

        # Fractional coefficients
        @test parse_formula("H1//2O") == Dict(:H => 1 // 2, :O => 1)
        @test parse_formula("Fe2//3O") == Dict(:Fe => 2 // 3, :O => 1)

        # Nested parentheses
        @test parse_formula("(NH4)2SO4") == Dict(:N => 2, :H => 8, :S => 1, :O => 4)

        # Charge extraction
        @test extract_charge("Ca+2") == 2
        @test extract_charge("SO4-2") == -2
        @test extract_charge("H2O") == 0
        @test extract_charge("H+") == 1
        @test extract_charge("OH-") == -1
    end

    @testset "Equation parsing" begin
        # Basic equation
        reactants, products, eq_sign = parse_equation("H2O = H+ + OH-")
        @test reactants == Dict("H2O" => 1)
        @test products == Dict("H+" => 1, "OH-" => 1)
        @test eq_sign == '='

        # With coefficients
        reactants, products, eq_sign = parse_equation("2H2O = 2H2 + O2")
        @test reactants == Dict("H2O" => 2)
        @test products == Dict("H2" => 2, "O2" => 1)

        # Forward arrow
        reactants, products, eq_sign = parse_equation("H2O → H+ + OH-")
        @test eq_sign == '→'

        # Empty left side (∅)
        reactants, products, eq_sign = parse_equation("∅ = H+ + OH-")
        @test isempty(reactants)
        @test length(products) == 2

        # parse_equation always returns 3-tuple
        result = parse_equation("CaCO3 = Ca+2 + CO3-2")
        @test length(result) == 3
        reac, prod, sign = result
        @test reac == Dict("CaCO3" => 1)
        @test prod == Dict("Ca+2" => 1, "CO3-2" => 1)
    end

    @testset "Formula formatting" begin
        # format_equation
        coeffs = Dict("H2O" => -2, "H2" => 2, "O2" => 1)
        eq = format_equation(coeffs)
        @test "2H2O = 2H2 + O2" == eq || "2H2O = O2 + 2H2" == eq

        # Charge balancing with electrons
        coeffs = Dict("Fe+2" => 2, "Fe+3" => 3)
        eq = format_equation(coeffs)
        @test occursin("e⁻", eq)

        # colored_equation does not throw
        @test_nowarn ChemistryLab.colored_equation("H2O = H+ + OH-")

        # add_parentheses_if_needed
        @test ChemistryLab.add_parentheses_if_needed("H2") == "H2"
        @test ChemistryLab.add_parentheses_if_needed("H+") == "(H+)"
    end

    @testset "Unicode conversion" begin
        @test phreeqc_to_unicode("Fe+3") == "Fe³⁺"
        @test phreeqc_to_unicode("SO4-2") == "SO₄²⁻"

        @test unicode_to_phreeqc("Fe³⁺") == "Fe+3"
        @test unicode_to_phreeqc("SO₄²⁻") == "SO4-2"

        # Fraction handling
        @test phreeqc_to_unicode("1//2H2O") == "½H₂O"
        @test unicode_to_phreeqc("½H₂O") == "1//2H2O"

        # Bidirectional round-trip
        for s in ["H2O", "Ca+2", "SO4-2", "OH-"]
            @test unicode_to_phreeqc(phreeqc_to_unicode(s)) == s
        end
    end

    @testset "Utility functions" begin
        # stoich_coef_round
        @test stoich_coef_round(2.0) == 2
        @test stoich_coef_round(1 / 2) == 1 // 2
        @test stoich_coef_round(0.3333333) ≈ 1 // 3

        # calculate_molar_mass
        atoms = Dict(:H => 2, :O => 1)
        @test calculate_molar_mass(atoms) ≈ 18.015u"g/mol"

        # to_mendeleev: oxide notation → element notation
        oxides = Dict(:C => 1, :S => 1)   # CaO + SiO2
        mendeleev = to_mendeleev(oxides)
        @test mendeleev[:Ca] == 1
        @test mendeleev[:O] == 3
        @test mendeleev[:Si] == 1

        # Multiple oxides
        oxides = Dict(:C => 2, :M => 1)   # 2CaO + MgO
        mendeleev = to_mendeleev(oxides)
        @test mendeleev[:Ca] == 2
        @test mendeleev[:Mg] == 1
        @test mendeleev[:O] == 3
    end
end
