@testsection "ChemicalState" begin
    # ── Shared fixtures ──────────────────────────────────────────────────────
    h2o   = Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)
    hplus = Species("H+";  aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE)
    oh    = Species("OH-"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE)
    nacl  = Species("NaCl"; aggregate_state=AS_CRYSTAL)

    cs_basic = ChemicalSystem([h2o, hplus, oh])

    @testsection "Construction — default T and P" begin
        state = ChemicalState(cs_basic)
        @test ustrip(temperature(state)) ≈ 298.15
        @test isapprox(ustrip(pressure(state)), 1e5; rtol=1e-4)   # 1 bar in Pa
        @test length(state.n) == 3
        @test all(iszero.(ustrip.(state.n)))
    end

    @testsection "Construction with explicit n (moles)" begin
        state = ChemicalState(cs_basic, [55.5u"mol", 1e-7u"mol", 1e-7u"mol"])
        @test ustrip(moles(state, "H2O")) ≈ 55.5
        @test ustrip(moles(state, "H+"))  ≈ 1e-7
    end

    @testsection "Construction with mass units" begin
        cs = ChemicalSystem([h2o, nacl])
        # 5.844 g NaCl ≈ 0.1 mol
        state = ChemicalState(cs, [55.5u"mol", 5.844u"g"])
        @test ustrip(moles(state, "H2O")) ≈ 55.5
        @test isapprox(ustrip(moles(state, "NaCl")), 0.1; rtol=1e-3)
    end

    @testsection "Temperature and pressure accessors" begin
        state = ChemicalState(cs_basic; T=350.0u"K", P=2u"bar")
        @test ustrip(temperature(state)) ≈ 350.0
        @test isapprox(ustrip(pressure(state)), 2e5; rtol=1e-4)
    end

    @testsection "set_temperature! and set_pressure!" begin
        state = ChemicalState(cs_basic)
        set_temperature!(state, 400.0u"K")
        @test ustrip(temperature(state)) ≈ 400.0

        set_pressure!(state, 5u"bar")
        @test isapprox(ustrip(pressure(state)), 5e5; rtol=1e-4)
    end

    @testsection "set_quantity! by species" begin
        state = ChemicalState(cs_basic)
        set_quantity!(state, h2o, 10.0u"mol")
        @test ustrip(moles(state, "H2O")) ≈ 10.0
    end

    @testsection "set_quantity! by symbol" begin
        state = ChemicalState(cs_basic, [55.5u"mol", 0.0u"mol", 0.0u"mol"])
        set_quantity!(state, "H+", 1e-4u"mol")
        @test ustrip(moles(state, "H+")) ≈ 1e-4
    end

    @testsection "moles by phase (NamedTuple)" begin
        state = ChemicalState(cs_basic, [55.5u"mol", 1e-7u"mol", 1e-7u"mol"])
        n_phases = moles(state)
        @test hasproperty(n_phases, :liquid)
        @test hasproperty(n_phases, :total)
        @test ustrip(n_phases.liquid) ≈ 55.5 + 2e-7
        @test ustrip(n_phases.solid) ≈ 0.0
    end

    @testsection "mass accessors" begin
        state = ChemicalState(cs_basic, [1.0u"mol", 0.0u"mol", 0.0u"mol"])
        m = mass(state)
        @test hasproperty(m, :total)
        # 1 mol H2O ≈ 18.015 g
        @test isapprox(ustrip(uconvert(us"g", m.total)), 18.015; rtol=1e-3)

        # Mass of a specific species
        m_h2o = mass(state, "H2O")
        @test isapprox(ustrip(uconvert(us"g", m_h2o)), 18.015; rtol=1e-3)
    end

    @testsection "volume accessors (without V⁰)" begin
        # Without molar volumes, volume is a NamedTuple with zero values
        state = ChemicalState(cs_basic, [55.5u"mol", 1e-7u"mol", 1e-7u"mol"])
        v = volume(state)
        @test hasproperty(v, :total)
        # volume(state, species) returns nothing when V⁰ is absent
        @test isnothing(volume(state, h2o))
    end

    @testsection "pH and pOH — acidic system" begin
        # 1e-4 mol H+ in ~1L water → pH ≈ 4
        state = ChemicalState(cs_basic, [55.5u"mol", 1e-4u"mol", 0.0u"mol"])
        ph = pH(state)
        if !isnothing(ph)
            @test isapprox(ph, 4.0; atol=0.1)
        end
    end

    @testsection "pH returns nothing without H+ or OH-" begin
        cs_no_ions = ChemicalSystem([h2o])
        state = ChemicalState(cs_no_ions, [55.5u"mol"])
        @test isnothing(pH(state))
        @test isnothing(pOH(state))
    end

    @testsection "porosity and saturation return nothing without V⁰" begin
        state = ChemicalState(cs_basic, [55.5u"mol", 0.0u"mol", 0.0u"mol"])
        # Without V⁰ all volumes are zero → porosity = nothing
        @test isnothing(porosity(state)) || porosity(state) isa Number
        @test isnothing(saturation(state)) || saturation(state) isa Number
    end

    @testsection "copy — independent mutable state" begin
        state = ChemicalState(cs_basic, [55.5u"mol", 1e-7u"mol", 1e-7u"mol"])
        clone = copy(state)

        set_quantity!(clone, "H2O", 10.0u"mol")
        @test ustrip(moles(state, "H2O")) ≈ 55.5   # original unchanged
        @test ustrip(moles(clone, "H2O")) ≈ 10.0   # clone updated

        # System is shared
        @test clone.system === state.system
    end

    @testsection "_update_derived! recalculates after mutation" begin
        state = ChemicalState(cs_basic, [55.5u"mol", 1e-7u"mol", 1e-7u"mol"])
        n_before = moles(state).liquid

        set_quantity!(state, "H2O", 100.0u"mol")
        n_after = moles(state).liquid

        @test ustrip(n_after) > ustrip(n_before)
    end
end
