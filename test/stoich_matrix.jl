using ChemistryLab
using Test
using OrderedCollections
using LinearAlgebra

@testsection "Stoichiometric Matrix" begin
    @testsection "Basic matrix construction" begin
        # Test simple species
        h2o = Species("H2O")
        hplus = Species("H+")
        oh = Species("OH-")
        species = [h2o, hplus, oh]
        
        # Test canonical matrix
        A, atoms = canonical_stoich_matrix(species; display=false)
        @test size(A) == (3, 3)  # H, O atoms and e- x 3 species
        @test atoms == [:H, :O, :Zz]
        @test A[1,1] == 2  # H2O has 2 H
        @test A[2,1] == 1  # H2O has 1 O
    end

    @testsection "Matrix with charged species" begin
        # Test with charged species
        species = [Species("Fe2+"), Species("Fe3+"), Species("e-")]
        A, indep_comp, dep_comp = stoich_matrix(species; display=false)
        
        # Check that charge is properly handled
        @test length(dep_comp) == 3  # Fe2+, Fe3+ and e-
        @test any(s -> charge(s) != 0, dep_comp)  # At least one species should be charged
    end

    @testsection "Reaction conversion" begin
        # Test conversion of stoichiometric matrix to reactions
        na2so4 = Species("Na2(SO4)")
        na = Species("Na+")
        so4 = Species("SO4-2")
        naso4 = Species("Na(SO4)-")
        species = [na2so4, na, so4, naso4]
        
        A, indep_comp, dep_comp = stoich_matrix(species; display=false)
        reactions = stoich_matrix_to_reactions(A, indep_comp, dep_comp; display=false)
        
        @test length(reactions) > 0
        @test all(r -> r isa Reaction, reactions)
    end

    @testsection "Mass-based calculations" begin
        # Test mass-based stoichiometric matrix
        h2o = Species("H2O")
        h2 = Species("H2")
        o2 = Species("O2")
        species = [h2o, h2, o2]
        
        A, atoms = canonical_stoich_matrix(species; display=false, mass=true)
        @test size(A) == (2, 3)
        
        # Mass conservation check (approximate due to floating point)
        h2o_idx = findfirst(s -> s == h2o, species)
        @test sum(A[:,h2o_idx]) ≈ 1.0 atol=1e-10
    end

    @testsection "CemSpecies handling" begin
        # Test handling of cement species
        species = [CemSpecies("C3S"), CemSpecies("CH"), CemSpecies("CSH")]
        A, atoms = canonical_stoich_matrix(species; display=false)
        @test size(A, 1) ≥ 2  # Should have at least Ca and Si components
        
        # Test that oxides are properly handled
        @test any(x -> x in [:C, :S, :H], atoms)  # Should contain basic oxide components
    end

    @testsection "Utility functions" begin
        # Test union_atoms
        d1 = OrderedDict(:Ca => 1, :O => 1)
        d2 = OrderedDict(:Si => 1, :O => 2)
        atoms = union_atoms([d1, d2])
        @test :Ca in atoms
        @test :Si in atoms
        @test :O in atoms
        
        # Test same_components functions
        species = [Species("CaCO3")]
        @test ChemistryLab.same_components(species) == atoms_charge
        
        cem_species = [CemSpecies("C3S")]
        @test ChemistryLab.same_components(cem_species) == oxides_charge
    end
end