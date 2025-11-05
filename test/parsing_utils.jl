using ChemistryLab
using Test

@testset "Parsing Utils" begin
    @testset "Character conversion" begin
        # Test superscript conversions
        @test ChemistryLab.super_to_normal("H₂O²⁺") == "H₂O2+"
        @test ChemistryLab.normal_to_super("2+") == "²⁺"
        
        # Test subscript conversions
        @test ChemistryLab.sub_to_normal("H₂O") == "H2O"
        @test ChemistryLab.normal_to_sub("H2O") == "H₂O"
        
        # Test all_normal_to_sub
        @test ChemistryLab.all_normal_to_sub("H2O") == "H₂O"
        @test ChemistryLab.all_normal_to_sub("SO4") == "SO₄"
    end

    @testset "Formula parsing" begin
        # Test basic formula parsing
        @test parse_formula("H2O") == Dict(:H => 2, :O => 1)
        @test parse_formula("Ca(OH)2") == Dict(:Ca => 1, :H => 2, :O => 2)
        
        # Test fractional coefficients
        @test parse_formula("H1//2O") == Dict(:H => 1//2, :O => 1)
        @test parse_formula("Fe2//3O") == Dict(:Fe => 2//3, :O => 1)
        
        # Test parentheses handling
        @test parse_formula("(NH4)2SO4") == Dict(:N => 2, :H => 8, :S => 1, :O => 4)
        
        # Test charge extraction
        @test extract_charge("Ca+2") == 2
        @test extract_charge("SO4-2") == -2
        @test extract_charge("H2O") == 0
    end

    @testset "Equation parsing" begin
        # Test basic equation parsing
        reactants, products, equal_sign = parse_equation("H2O = H+ + OH-")
        @test reactants == Dict("H2O" => 1)
        @test products == Dict("H+" => 1, "OH-" => 1)
        @test equal_sign == '='
        
        # Test with coefficients
        reactants, products, equal_sign = parse_equation("2H2O = 2H2 + O2")
        @test reactants == Dict("H2O" => 2)
        @test products == Dict("H2" => 2, "O2" => 1)
        
        # Test with different reaction arrows
        reactants, products, equal_sign = parse_equation("H2O → H+ + OH-")
        @test equal_sign == '→'
        
        # Test empty sides
        reactants, products, equal_sign = parse_equation("∅ = H+ + OH-")
        @test isempty(reactants)
        @test length(products) == 2
    end

    @testset "Formula formatting" begin
        # Test colored formula
        # @test occursin("H₂O", colored_formula("H2O"))
        # @test occursin("Ca²⁺", colored_formula("Ca+2"))
        
        # Test equation formatting
        coeffs = Dict("H2O" => -2, "H2" => 2, "O2" => 1)
        eq = format_equation(coeffs)
        @test "2H2O = 2H2 + O2" == eq || "2H2O = O2 + 2H2" == eq
        
        # Test with charge balancing
        coeffs = Dict("Fe+2" => 2, "Fe+3" => 3)
        eq = format_equation(coeffs)
        @test occursin("e⁻", eq)  # Should add electron to balance charge
    end

    @testset "Unicode conversion" begin
        # Test phreeqc to unicode conversion
        @test phreeqc_to_unicode("Fe+3") == "Fe³⁺"
        @test phreeqc_to_unicode("SO4-2") == "SO₄²⁻"
        
        # Test unicode to phreeqc conversion
        @test unicode_to_phreeqc("Fe³⁺") == "Fe+3"
        @test unicode_to_phreeqc("SO₄²⁻") == "SO4-2"
        
        # Test fraction handling
        @test phreeqc_to_unicode("1//2H2O") == "½H₂O"
        @test unicode_to_phreeqc("½H₂O") == "1//2H2O"
    end

    @testset "Utility functions" begin
        # Test stoich_coef_round
        @test stoich_coef_round(2.0) == 2
        @test stoich_coef_round(1/2) == 1//2
        @test stoich_coef_round(0.3333333) ≈ 1//3
        
        # Test calculate_molar_mass
        atoms = Dict(:H => 2, :O => 1)
        @test calculate_molar_mass(atoms) ≈ 18.015u"g/mol"
        
        # # Test merge_sum_dicts
        # d1 = Dict(:H => 2, :O => 1)
        # d2 = Dict(:H => 1, :O => 2)
        # merged = ChemistryLab.merge_sum_dicts([d1, d2])
        # @test merged[:H] == 3
        # @test merged[:O] == 3
    end

    @testset "Cement to Mendeleev conversion" begin
        # Test oxide conversion
        oxides = Dict(:C => 1, :S => 1)  # CaO and SiO2
        mendeleev = to_mendeleev(oxides)
        @test mendeleev[:Ca] == 1
        @test mendeleev[:O] == 3
        @test mendeleev[:Si] == 1
        
        # Test multiple oxides
        oxides = Dict(:C => 2, :M => 1)  # 2CaO + MgO
        mendeleev = to_mendeleev(oxides)
        @test mendeleev[:Ca] == 2
        @test mendeleev[:Mg] == 1
        @test mendeleev[:O] == 3
    end
end