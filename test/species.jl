@testsection "Species" begin
    # Test basic Species construction
    @test_nowarn Species("H2O")
    water = Species("H2O", name = "Water", symbol = "H2O", aggregate_state = AS_AQUEOUS)
    @test name(water) == "Water"
    @test symbol(water) == "H2O"
    @test aggregate_state(water) == AS_AQUEOUS

    # Test atoms and composition
    @test atoms(water) == Dict(:H => 2, :O => 1)
    @test charge(water) == 0

    # Test properties
    @test water.M ≈ 18.015u"g/mol"  # Molar mass of water

    # Test property manipulation
    water[:custom_prop] = 42
    @test water[:custom_prop] == 42
    @test haskey(water, :custom_prop)
    @test !haskey(water, :nonexistent)

    # Test species equality (based on formula + aggregate_state + class)
    water2 = Species("H2O", name = "Water2", aggregate_state = AS_AQUEOUS)
    @test water == water2

    # Test species inequality: same formula, different aggregate_state
    vapour = Species("H₂O"; name = "Vapour", symbol = "H₂O⤴", aggregate_state = AS_GAS, class = SC_GASFLUID)
    @test vapour != water

    # Test ionic species
    nacl = Species("Na+")
    @test charge(nacl) == 1
    @test atoms(nacl) == Dict(:Na => 1)

    # Test complex formula
    calcium_carbonate = Species("CaCO3", aggregate_state = AS_CRYSTAL)
    @test atoms(calcium_carbonate) == Dict(:Ca => 1, :C => 1, :O => 3)
    @test aggregate_state(calcium_carbonate) == AS_CRYSTAL

    # Test aqueous solvent class
    h2o_solvent = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    @test class(h2o_solvent) == SC_AQSOLVENT

    # Test aqueous solute class
    na_plus = Species("Na+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    @test class(na_plus) == SC_AQSOLUTE

    # Test component class
    sio2 = Species("SiO2"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    @test class(sio2) == SC_COMPONENT

    # Test mendeleev_filter (not exported — use ChemistryLab.mendeleev_filter)
    valid_species = Species("CaCO3"; aggregate_state = AS_CRYSTAL)
    @test !isnothing(ChemistryLab.mendeleev_filter(valid_species))
    invalid_species = Species("Xx"; aggregate_state = AS_UNDEF)
    @test isnothing(ChemistryLab.mendeleev_filter(invalid_species))

    # Test apply: double all stoichiometric coefficients
    doubled = apply(x -> x * 2, water)
    @test atoms(doubled) == Dict(:H => 4, :O => 2)

    # Test aqueous dissolved species with @ suffix
    @test_nowarn Species("CO2@"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    co2_aq = Species("CO2@"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    @test aggregate_state(co2_aq) == AS_AQUEOUS

    # Test CemSpecies construction
    c3s = CemSpecies("C3S")
    @test !isempty(atoms(c3s))
    @test :Ca in keys(atoms(c3s))
    @test :Si in keys(atoms(c3s))

    # CemSpecies oxide composition
    ch = CemSpecies("CH")
    ch_oxides = oxides(ch)
    @test :C in keys(ch_oxides)   # CaO
    @test :H in keys(ch_oxides)   # H2O

    # CemSpecies components returns oxides
    @test oxides(c3s) == components(c3s)

    # CemSpecies hash consistency
    c3s_copy = CemSpecies("C3S")
    @test c3s == c3s_copy
    @test hash(c3s) == hash(c3s_copy)
end

@testsection "with_class" begin
    # Base case: SC_COMPONENT → SC_SSENDMEMBER
    s = Species("CaCO3"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    s2 = with_class(s, SC_SSENDMEMBER)

    @test class(s2) == SC_SSENDMEMBER
    @test class(s) == SC_COMPONENT          # original unchanged

    # All other fields preserved
    @test name(s2) == name(s)
    @test symbol(s2) == symbol(s)
    @test formula(s2) == formula(s)
    @test aggregate_state(s2) == aggregate_state(s)
    @test properties(s2) === properties(s)  # same dict object (shared by reference)

    # Round-trip: any class can be set
    s3 = with_class(s2, SC_AQSOLUTE)
    @test class(s3) == SC_AQSOLUTE
    @test class(s2) == SC_SSENDMEMBER       # s2 still unchanged

    # Works with a species carrying properties
    sp = Species("Ca+2"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    sp[:custom] = 42
    sp2 = with_class(sp, SC_UNDEF)
    @test class(sp2) == SC_UNDEF
    @test sp2[:custom] == 42

    # The result can be used in a SolidSolutionPhase
    em1 = with_class(
        Species("Em1"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT),
        SC_SSENDMEMBER,
    )
    em2 = with_class(
        Species("Em2"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT),
        SC_SSENDMEMBER,
    )
    @test_nowarn SolidSolutionPhase("SS", [em1, em2])
end
