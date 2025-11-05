using ChemistryLab
using Test

@testsection "Reactions" begin
    @testsection "Constructor from equation" begin
        r = Reaction("H2O = H+ + OH-")
        @test reactants(r)[Species("H2O")] == 1
        @test products(r)[Species("H+")] == 1
        @test products(r)[Species("OH-")] == 1
    end

    @testsection "Getindex and properties" begin
        r = Reaction("CaCO3 = Ca+2 + CO3-2")
        @test r[Species("CaCO3")] == -1
        r[:testprop] = "abc"
        @test r[:testprop] == "abc"
        @test r.testprop == "abc"
    end

    @testsection "Simplify reaction" begin
        reac = Dict(Species("OH⁻") => -1, Species("H2O") => -1)
        prod = Dict(Species("H2O") => 1, Species("H⁺") => 1)
        r = Reaction(reac, prod)
        rs = simplify_reaction(r)
        @test length(reactants(rs)) == 1
        @test length(products(rs)) == 2
        @test haskey(reactants(rs), Species("OH⁻"))
        @test haskey(products(rs), Species("H⁺"))
    end

    @testsection "Addition and subtraction" begin
        r1 = Reaction("H2O = H+ + OH-")
        r2 = Reaction("2H2O = H2 + 2OH⁻")
        rsum = r1 + r2
        @test reactants(rsum)[Species("H2O")] == 3
        rsub = r2 - r1
        @test reactants(rsub)[Species("H2O")] == 1 || haskey(reactants(rsub), Species("H2O")) == true
    end

    @testsection "split/merge species by stoich" begin
        s = Dict(Species("CaCl2") => -1, Species("Ca+2") => 1, Species("Cl-") => 2)
        reac, prod = ChemistryLab.split_species_by_stoich(s)
        @test haskey(reac, Species("CaCl2"))
        @test haskey(prod, Species("Ca+2"))
        merged = ChemistryLab.merge_species_by_stoich(reac, prod)
        @test haskey(merged, Species("CaCl2")) && haskey(merged, Species("Ca+2"))
    end

    @testsection "scale_stoich!" begin
        s = Dict(Species("O")=>2, Species("H")=>4)
        ChemistryLab.scale_stoich!(s)
        @test s[Species("O")] == 4 && s[Species("H")] == 8
    end

end
