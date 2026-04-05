@testsection "Thermodynamics" begin
    @testsection "THERMO_MODELS registry" begin
        # Known models must be present
        @test haskey(THERMO_MODELS, :cp_ft_equation)
        @test haskey(THERMO_MODELS, :logk_fpt_function)

        # cp_ft_equation must have Cp, H, S, G expressions
        model = THERMO_MODELS[:cp_ft_equation]
        @test haskey(model, :Cp)
        @test haskey(model, :H)
        @test haskey(model, :S)
        @test haskey(model, :G)
        @test haskey(model, :units)

        # logk_fpt_function must have logKr
        model_logk = THERMO_MODELS[:logk_fpt_function]
        @test haskey(model_logk, :logKr)
    end

    @testsection "THERMO_FACTORIES populated at __init__" begin
        @test !isempty(THERMO_FACTORIES)
        @test haskey(THERMO_FACTORIES, :cp_ft_equation)
        @test haskey(THERMO_FACTORIES, :logk_fpt_function)

        factories = THERMO_FACTORIES[:cp_ft_equation]
        @test haskey(factories, :Cp)
        @test haskey(factories, :H)
        @test haskey(factories, :S)
        @test factories[:Cp] isa ThermoFactory
    end

    @testsection "SymbolicFunc from constant Number" begin
        tf = SymbolicFunc(42.0)
        @test tf isa SymbolicFunc
        @test tf() ≈ 42.0
    end

    @testsection "SymbolicFunc from constant Quantity" begin
        tf = SymbolicFunc(100.0u"J/mol")
        @test tf isa SymbolicFunc
        val = tf(; unit = true)
        @test val isa AbstractQuantity
    end

    @testsection "SymbolicFunc from expression" begin
        # Simple linear function: f(T) = 2*T
        tf = SymbolicFunc(:(2 * T), [:T])
        @test tf isa SymbolicFunc
        @test tf(; T = 3.0) ≈ 6.0
        @test tf(; T = 10.0) ≈ 20.0
    end

    @testsection "SymbolicFunc arithmetic" begin
        tf1 = SymbolicFunc(:(2 * T), [:T])
        tf2 = SymbolicFunc(:(3 * T), [:T])

        # Addition
        tf_sum = tf1 + tf2
        @test tf_sum(; T = 1.0) ≈ 5.0

        # Subtraction
        tf_diff = tf2 - tf1
        @test tf_diff(; T = 1.0) ≈ 1.0

        # Scalar multiplication
        tf_scaled = tf1 * 2.0
        @test tf_scaled(; T = 1.0) ≈ 4.0

        # Unary negation
        tf_neg = -tf1
        @test tf_neg(; T = 1.0) ≈ -2.0
    end

    @testsection "ThermoFactory construction and call" begin
        factory = ThermoFactory(:(a * T + b), [:T])
        @test factory isa ThermoFactory
        @test haskey(factory.params, :a)
        @test haskey(factory.params, :b)
        @test haskey(factory.vars, :T)

        # Instantiate with parameters
        tf = factory(; a = 2.0, b = 5.0)
        @test tf isa SymbolicFunc
        @test tf(; T = 10.0) ≈ 25.0   # 2*10 + 5
    end

    @testsection "add_thermo_model" begin
        model_name = :test_linear_model
        Cpexpr = :(c₀ + c₁ * T)

        # Register a new model from a Cp expression
        @test_nowarn add_thermo_model(model_name, Cpexpr)

        @test haskey(THERMO_MODELS, model_name)
        @test haskey(THERMO_FACTORIES, model_name)

        factories = THERMO_FACTORIES[model_name]
        @test haskey(factories, :Cp)
        @test haskey(factories, :H)   # integrated automatically

        # Clean up: remove the test model so it doesn't pollute other tests
        delete!(THERMO_MODELS, model_name)
        delete!(THERMO_FACTORIES, model_name)
    end

    @testsection "build_thermo_functions" begin
        # Build thermodynamic functions from cp_ft_equation with minimal parameters
        params = [
            :a₀ => 75.0u"J/(mol*K)",
            :a₁ => 0.0u"J/(mol*K^2)",
            :a₂ => 0.0u"J*K/mol",
            :a₃ => 0.0u"J/(mol*K^0.5)",
            :a₄ => 0.0u"J/(mol*K^3)",
            :a₅ => 0.0u"J/(mol*K^4)",
            :a₆ => 0.0u"J/(mol*K^5)",
            :a₇ => 0.0u"J*K^2/mol",
            :a₈ => 0.0u"J/mol",
            :a₉ => 0.0u"J/(mol*K^1.5)",
            :a₁₀ => 0.0u"J/(mol*K)",
            :T => 298.15u"K",
            :S⁰ => 70.0u"J/(mol*K)",
            :ΔfH⁰ => -285.8u"kJ/mol",
            :ΔfG⁰ => -237.1u"kJ/mol",
        ]
        thermo = build_thermo_functions(:cp_ft_equation, params)

        @test haskey(thermo, :Cp⁰)
        @test haskey(thermo, :ΔₐH⁰)
        @test haskey(thermo, :S⁰)
        @test haskey(thermo, :ΔₐG⁰)

        @test thermo[:Cp⁰] isa SymbolicFunc

        # Cp at Tref should be ≈ 75 J/(mol·K) (only constant term a₀)
        cp_val = thermo[:Cp⁰](; T = 298.15)
        @test isapprox(cp_val, 75.0; rtol = 1.0e-6)
    end
end
