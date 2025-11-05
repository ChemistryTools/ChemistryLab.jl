@testsection "Formulas" begin
    # Basic parsing and composition
    f = Formula("H2O")
    @test composition(f) == Dict(:H => 2, :O => 1)
    @test charge(f) == 0
    @test f[:H] == 2
    @test f[:He] == 0 # missing element returns zero
    @test length(f) == 2
    @test f == Formula("H2O")

    # Charge handling
    na_plus = Formula("Na+")
    @test charge(na_plus) == 1
    @test composition(na_plus) == Dict(:Na => 1)

    # Build from composition and ensure expr/composition consistency
    comp = Dict(:Ca => 1, :C => 1, :O => 3)
    fc = Formula(comp)
    @test composition(fc) == comp
    @test fc[:Ca] == 1

    # Arithmetic on formulas
    f2 = f * 2
    @test composition(f2) == Dict(:H => 4, :O => 2)
    f_half = f / 2
    @test composition(f_half) == Dict(:H => 1, :O => 0.5)

    # Addition with AtomGroup and Symbol
    f_plus_na = f + AtomGroup(:Na)
    @test composition(f_plus_na)[:Na] == 1
    f_plus_ag = f + AtomGroup(2, :Cl)
    @test composition(f_plus_ag)[:Cl] == 2

    # AtomGroup addition
    agsum = AtomGroup(:H) + AtomGroup(2, :H)
    @test isa(agsum, Formula)
    @test composition(agsum)[:H] == 3

    # convert: ensure numeric types change
    f_int = Formula(Dict(:H => 2, :O => 1))
    f_float = convert(Float64, f_int)
    vals = collect(values(composition(f_float)))
    @test all(x->x isa Float64, vals)

    # apply: apply an identity function should preserve composition values
    f_applied = apply(x->x, f_int)
    @test composition(f_applied) == composition(f_int)

    # check_mendeleev should throw on invalid atom
    @test_throws ErrorException check_mendeleev(Formula("Xx"))

end
