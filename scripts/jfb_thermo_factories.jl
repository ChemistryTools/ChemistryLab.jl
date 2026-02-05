using Revise
using ChemistryLab
using DynamicQuantities
using ModelingToolkit
using SymbolicNumericIntegration

# Construction
expr = :(α + β * T + γ * log(T))
factory = ThermoFactory(expr, [:T])
params = (α=210.0, β=0.0, γ=-3.07e6, T=298.15)
f = factory(; params...)
# or directly
f = factory(α=210.0, β=0.0, γ=-3.07e6, T=298.15)
f(T = 300)
f()
f(T = 298.15)

params = [:a₀ => 210.0, :a₁ => 0.12, :a₂ => -3.07e6, :a₃ => 0.0, :T => 298.15, :Cp⁰ => 210.0, :ΔfH⁰ => -2723484.33, :S⁰ => 140, :ΔfG⁰ => -2480808.197, :V⁰ => 7.84]

# For manual testing but the best method consists in completing models thanks to `add_thermo_model`
factoryCp = THERMO_FACTORIES[:cp_ft_equation][:Cp]
Cp⁰ = factoryCp(; params...)
H = integrate(Cp⁰.symbolic, only(@variables T); symbolic = true, detailed = false)
fH = ThermoFunction(H; T=298.15)

dtf = build_thermo_functions(:cp_ft_equation, params)
for (k, v) in dtf
    println(k)
    display(v)
    println()
end
