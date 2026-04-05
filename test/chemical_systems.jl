@testsection "ChemicalSystem" begin
    # ── Shared fixtures ──────────────────────────────────────────────────────
    h2o = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    hplus = Species("H+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    oh = Species("OH-"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    nacl = Species("NaCl"; aggregate_state = AS_CRYSTAL)
    co2g = Species("CO2"; aggregate_state = AS_GAS, class = SC_GASFLUID)
    sio2 = Species("SiO2"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)

    @testsection "Basic construction" begin
        cs = ChemicalSystem([h2o, hplus, oh])
        @test length(cs) == 3
        @test size(cs) == (3,)
    end

    @testsection "Indexing" begin
        cs = ChemicalSystem([h2o, hplus, oh])
        @test cs[1] == h2o
        @test cs["H2O"] == h2o
        @test cs["H+"] == hplus
        @test_throws KeyError cs["nonexistent"]
    end

    @testsection "Iteration" begin
        cs = ChemicalSystem([h2o, hplus, oh])
        collected = collect(cs)
        @test length(collected) == 3
        @test h2o in collected
    end

    @testsection "Views by aggregate state" begin
        cs = ChemicalSystem([h2o, hplus, nacl, co2g])
        @test length(aqueous(cs)) == 2   # h2o + hplus
        @test length(crystal(cs)) == 1   # nacl
        @test length(gas(cs)) == 1       # co2g
        @test all(s -> aggregate_state(s) == AS_AQUEOUS, aqueous(cs))
        @test all(s -> aggregate_state(s) == AS_CRYSTAL, crystal(cs))
        @test all(s -> aggregate_state(s) == AS_GAS, gas(cs))
    end

    @testsection "Views by class" begin
        cs = ChemicalSystem([h2o, hplus, oh, co2g, sio2])
        @test length(solutes(cs)) == 2      # hplus + oh
        @test class(solvent(cs)) == SC_AQSOLVENT
        @test length(gasfluid(cs)) == 1     # co2g
        @test length(components(cs)) == 1   # sio2
    end

    @testsection "Construction with explicit primaries" begin
        sp = [h2o, hplus, oh]
        cs = ChemicalSystem(sp, [h2o])
        @test length(cs.SM.primaries) == 1
        @test symbol(cs.SM.primaries[1]) == "H2O"
    end

    @testsection "Construction from string primaries" begin
        sp = [h2o, hplus, oh]
        cs = ChemicalSystem(sp, ["H2O", "H+"])
        primaries_syms = symbol.(cs.SM.primaries)
        @test "H2O" in primaries_syms
        @test "H+" in primaries_syms
        @test "OH-" ∉ primaries_syms
    end

    @testsection "Pre-computed matrices" begin
        cs = ChemicalSystem([h2o, hplus, oh])
        @test !isnothing(cs.CSM)   # canonical matrix built
        @test !isnothing(cs.SM)    # primary-relative matrix built
        @test size(cs.CSM.A, 2) == 3   # 3 species → 3 columns
    end

    @testsection "Reactions and get_reaction" begin
        r_diss = Reaction("H2O = H+ + OH-"; symbol = "water_diss")
        cs = ChemicalSystem([h2o, hplus, oh]; reactions = [r_diss])
        @test length(cs.reactions) == 1
        @test symbol(get_reaction(cs, "water_diss")) == "water_diss"
        @test_throws KeyError get_reaction(cs, "nonexistent_rxn")
    end

    @testsection "merge two systems" begin
        cs1 = ChemicalSystem(
            [h2o, hplus];
            reactions = [Reaction("H2O = H+ + OH-"; symbol = "water_diss")],
        )
        cs2 = ChemicalSystem(
            [oh, h2o];   # h2o duplicate — cs1 wins
            reactions = [Reaction("CO2 + H2O = H2CO3"; symbol = "co2_hyd")],
        )
        cs = merge(cs1, cs2)
        @test length(cs) == 3             # H2O, H+, OH- (no duplicate H2O)
        @test length(cs.reactions) == 2   # both reactions present
        @test cs["H2O"] === cs1["H2O"]    # cs1 wins on conflict
    end

    @testsection "merge variadic" begin
        cs1 = ChemicalSystem([h2o])
        cs2 = ChemicalSystem([hplus])
        cs3 = ChemicalSystem([oh])
        cs = merge(cs1, cs2, cs3)
        @test length(cs) == 3
    end
end
