@testset "Species" begin
    # Test basic Species construction
    @test_nowarn Species("H2O")
    water = Species("H2O", name="Water", symbol="H2O", aggregate_state=AS_AQUEOUS)
    @test name(water) == "Water"
    @test symbol(water) == "H2O"
    @test aggregate_state(water) == AS_AQUEOUS
    
    # Test atoms and composition
    @test atoms(water) == OrderedDict(:H => 2, :O => 1)
    @test charge(water) == 0
    
    # Test properties
    @test ustrip(water.M) ≈ 18.015 atol=0.001  # Molar mass of water
    
    # Test property manipulation
    water[:custom_prop] = 42
    @test water[:custom_prop] == 42
    @test haskey(water, :custom_prop)
    
    # Test species equality
    water2 = Species("H2O", name="Water2", aggregate_state=AS_AQUEOUS)
    @test water == water2  # Should be equal because formula and state are the same

    # Test species difference
    fH₂O = 2*:H + :O
    vapour = Species(fH₂O; name="Vapour", symbol="H₂O⤴", aggregate_state=AS_GAS, class=SC_GASFLUID)
    vapour != water # false since aggregate_state or class are different despite atoms are identical

    
    # Test ionic species
    nacl = Species("Na+")
    @test charge(nacl) == 1
    @test atoms(nacl) == OrderedDict(:Na => 1)
    
    # Test complex formula
    calcium_carbonate = Species("CaCO3", aggregate_state=AS_CRYSTAL)
    @test atoms(calcium_carbonate) == OrderedDict(:Ca => 1, :C => 1, :O => 3)
    @test aggregate_state(calcium_carbonate) == AS_CRYSTAL
end
#IUPAC TODO
