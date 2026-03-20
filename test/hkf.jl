using ForwardDiff

@testsection "HKF water properties (HGK EOS)" begin
    wtp = water_thermo_props(298.15, 1.0e5)

    # Liquid water density at 25°C, 1 bar ≈ 997 kg/m³
    @test isapprox(wtp.D, 997.0; rtol = 1e-2)

    # DT should be negative (water expands when heated)
    @test wtp.DT < 0

    # DP should be positive (density increases under pressure)
    @test wtp.DP > 0
end

@testsection "HKF water electrostatics (Johnson-Norton)" begin
    wtp = water_thermo_props(298.15, 1.0e5)
    wep = water_electro_props_jn(298.15, 1.0e5, wtp)

    # Dielectric constant of water at 25°C ≈ 78.4
    @test isapprox(wep.epsilon, 78.4; rtol = 1e-2)

    # Born function Z = -1/ε
    @test isapprox(wep.bornZ, -1 / wep.epsilon; rtol = 1e-6)

    # Y = εT/ε² should be negative (ε decreases with T)
    @test wep.bornY < 0
end

@testsection "HKF g-function (Shock 1992)" begin
    wtp = water_thermo_props(298.15, 1.0e5)
    gs  = hkf_g_function(298.15, 1.0e5, wtp)

    # g ≈ 0 at reference conditions (density ~997, well inside validity range)
    @test isapprox(gs.g, 0.0; atol = 1e-2)

    # Outside validity range → zero state
    wtp_bad = WaterThermoProps(1100.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    gs_bad  = hkf_g_function(298.15, 1.0e5, wtp_bad)
    @test iszero(gs_bad.g)
end

@testsection "HKF thermo functions — Al(OH)2+" begin
    # Parameters for Al(OH)2+ from aq17-thermofun.json, converted to SI
    # SUPCRT values × HKF_SI_CONVERSIONS factors
    params = [
        :a1   =>  0.24940000474453 * 4.184e-5,
        :a2   => -169.08999633789  * 4.184,
        :a3   =>  6.4145998954773  * 4.184e-5,
        :a4   => -27091.0          * 4.184,
        :c1   =>  16.743900299072  * 4.184,
        :c2   => -10465.0          * 4.184,
        :wref =>  53240.0          * 4.184,
        :z    =>  1.0,
        :ΔₐG⁰ => -898292.0u"J/mol",
        :ΔₐH⁰ => -995581.0u"J/mol",
        :S⁰   => -27.53u"J/(mol*K)",
        :T    => 298.15u"K",
        :P    => 1.0e5u"Pa",
    ]

    thermo = build_thermo_functions(:solute_hkf88_reaktoro, params)

    @test haskey(thermo, :Cp⁰)
    @test haskey(thermo, :ΔₐH⁰)
    @test haskey(thermo, :S⁰)
    @test haskey(thermo, :ΔₐG⁰)
    @test haskey(thermo, :V⁰)

    @test thermo[:Cp⁰]  isa NumericFunc
    @test thermo[:ΔₐG⁰] isa NumericFunc

    # vars must be (:T, :P)
    @test thermo[:Cp⁰].vars == (:T, :P)

    # refs must store T and P as Quantities in SI
    @test haskey(thermo[:Cp⁰].refs, :T)
    @test haskey(thermo[:Cp⁰].refs, :P)
    @test thermo[:Cp⁰].refs.T isa AbstractQuantity
    @test thermo[:Cp⁰].refs.P isa AbstractQuantity
    @test ustrip(thermo[:Cp⁰].refs.T) ≈ 298.15
    @test ustrip(thermo[:Cp⁰].refs.P) ≈ 1.0e5

    # Cp at (Tr, Pr) must match database value ≈ 40.87 J/(mol·K)
    @test isapprox(thermo[:Cp⁰](; T = 298.15, P = 1.0e5), 40.87; rtol = 1e-2)

    # Calling without kwargs must use refs and give the same result
    @test isapprox(thermo[:Cp⁰](), thermo[:Cp⁰](; T = 298.15, P = 1.0e5); rtol = 1e-10)

    # G at (Tr, Pr) must match Gf
    @test isapprox(thermo[:ΔₐG⁰](; T = 298.15, P = 1.0e5), -898292.0; rtol = 1e-3)
    @test isapprox(thermo[:ΔₐG⁰](), -898292.0; rtol = 1e-3)

    # unit=true returns a Quantity
    @test thermo[:Cp⁰](; T = 298.15, P = 1.0e5, unit = true) isa AbstractQuantity
end

@testsection "HKF ForwardDiff compatibility" begin
    params = [
        :a1   =>  0.24940000474453 * 4.184e-5,
        :a2   => -169.08999633789  * 4.184,
        :a3   =>  6.4145998954773  * 4.184e-5,
        :a4   => -27091.0          * 4.184,
        :c1   =>  16.743900299072  * 4.184,
        :c2   => -10465.0          * 4.184,
        :wref =>  53240.0          * 4.184,
        :z    =>  1.0,
        :ΔₐG⁰ => -898292.0u"J/mol",
        :ΔₐH⁰ => -995581.0u"J/mol",
        :S⁰   => -27.53u"J/(mol*K)",
        :T    => 298.15u"K",
        :P    => 1.0e5u"Pa",
    ]
    thermo = build_thermo_functions(:solute_hkf88_reaktoro, params)

    # ∂G/∂T = -S (thermodynamic identity)
    dG_dT = ForwardDiff.derivative(T -> thermo[:ΔₐG⁰](; T = T, P = 1.0e5), 298.15)
    S_val  = thermo[:S⁰](; T = 298.15, P = 1.0e5)
    @test isapprox(dG_dT, -S_val; rtol = 1e-3)
end

@testsection "NumericFunc arithmetic — scalar" begin
    params = [
        :a1   =>  0.24940000474453 * 4.184e-5,
        :a2   => -169.08999633789  * 4.184,
        :a3   =>  6.4145998954773  * 4.184e-5,
        :a4   => -27091.0          * 4.184,
        :c1   =>  16.743900299072  * 4.184,
        :c2   => -10465.0          * 4.184,
        :wref =>  53240.0          * 4.184,
        :z    =>  1.0,
        :ΔₐG⁰ => -898292.0u"J/mol",
        :ΔₐH⁰ => -995581.0u"J/mol",
        :S⁰   => -27.53u"J/(mol*K)",
        :T    => 298.15u"K",
        :P    => 1.0e5u"Pa",
    ]
    thermo = build_thermo_functions(:solute_hkf88_reaktoro, params)
    G = thermo[:ΔₐG⁰]
    T0, P0 = 298.15, 1.0e5
    G0 = G(; T = T0, P = P0)

    # Unary negation
    neg = -G
    @test neg isa NumericFunc
    @test neg(; T = T0, P = P0) ≈ -G0

    # Scalar multiplication (stoichiometric coefficient)
    scaled = 2 * G
    @test scaled isa NumericFunc
    @test scaled(; T = T0, P = P0) ≈ 2 * G0

    scaled2 = G * 3.0
    @test scaled2 isa NumericFunc
    @test scaled2(; T = T0, P = P0) ≈ 3.0 * G0

    # Scalar division
    div_hkf = G / 2.0
    @test div_hkf isa NumericFunc
    @test div_hkf(; T = T0, P = P0) ≈ G0 / 2.0

    # Scalar addition/subtraction (offset in same unit)
    offset = 1000.0
    add_hkf = G + offset
    @test add_hkf(; T = T0, P = P0) ≈ G0 + offset
    sub_hkf = G - offset
    @test sub_hkf(; T = T0, P = P0) ≈ G0 - offset
end

@testsection "NumericFunc arithmetic — HKF op HKF" begin
    params = [
        :a1   =>  0.24940000474453 * 4.184e-5,
        :a2   => -169.08999633789  * 4.184,
        :a3   =>  6.4145998954773  * 4.184e-5,
        :a4   => -27091.0          * 4.184,
        :c1   =>  16.743900299072  * 4.184,
        :c2   => -10465.0          * 4.184,
        :wref =>  53240.0          * 4.184,
        :z    =>  1.0,
        :ΔₐG⁰ => -898292.0u"J/mol",
        :ΔₐH⁰ => -995581.0u"J/mol",
        :S⁰   => -27.53u"J/(mol*K)",
        :T    => 298.15u"K",
        :P    => 1.0e5u"Pa",
    ]
    thermo = build_thermo_functions(:solute_hkf88_reaktoro, params)
    G  = thermo[:ΔₐG⁰]
    H  = thermo[:ΔₐH⁰]
    T0, P0 = 298.15, 1.0e5
    G0 = G(; T = T0, P = P0)
    H0 = H(; T = T0, P = P0)

    sum_hkf = G + H
    @test sum_hkf isa NumericFunc
    @test sum_hkf(; T = T0, P = P0) ≈ G0 + H0

    diff_hkf = G - H
    @test diff_hkf isa NumericFunc
    @test diff_hkf(; T = T0, P = P0) ≈ G0 - H0

    # Linear combination (reaction-like)
    lc = 2 * G + (-1) * H
    @test lc isa NumericFunc
    @test lc(; T = T0, P = P0) ≈ 2 * G0 - H0
end

@testsection "NumericFunc arithmetic — cross-type with SymbolicFunc" begin
    params = [
        :a1   =>  0.24940000474453 * 4.184e-5,
        :a2   => -169.08999633789  * 4.184,
        :a3   =>  6.4145998954773  * 4.184e-5,
        :a4   => -27091.0          * 4.184,
        :c1   =>  16.743900299072  * 4.184,
        :c2   => -10465.0          * 4.184,
        :wref =>  53240.0          * 4.184,
        :z    =>  1.0,
        :ΔₐG⁰ => -898292.0u"J/mol",
        :ΔₐH⁰ => -995581.0u"J/mol",
        :S⁰   => -27.53u"J/(mol*K)",
        :T    => 298.15u"K",
        :P    => 1.0e5u"Pa",
    ]
    thermo = build_thermo_functions(:solute_hkf88_reaktoro, params)
    G  = thermo[:ΔₐG⁰]
    S  = thermo[:S⁰]
    T0, P0 = 298.15, 1.0e5
    G0 = G(; T = T0, P = P0)
    S0 = S(; T = T0, P = P0)

    # SymbolicFunc representing a constant (must carry same unit as G for +/-)
    tf_const = SymbolicFunc(-898292.0u"J/mol")

    add_cross = G + tf_const
    @test add_cross isa NumericFunc
    @test add_cross(; T = T0, P = P0) ≈ G0 + (-898292.0)

    diff_cross = tf_const - G
    @test diff_cross isa NumericFunc
    @test diff_cross(; T = T0, P = P0) ≈ -898292.0 - G0

    # Division of HKF by a SymbolicFunc of T (as in logK computation)
    R_log10 = ustrip(DynamicQuantities.Constants.R) * log(10)
    tf_T = R_log10 * SymbolicFunc(:T; units = [:T => "K"])
    logK = -G / tf_T
    @test logK isa NumericFunc
    # logK = -G / (R*ln10*T), numerical check
    expected_logK = -G0 / (R_log10 * T0)
    @test isapprox(logK(; T = T0, P = P0), expected_logK; rtol = 1e-6)
end
