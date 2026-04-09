using PrecompileTools: @compile_workload

# Exercise the ThermoFactory → SymbolicFunc compilation pipeline at precompile time
# to reduce first-call latency at runtime.
#
# IMPORTANT: uparse(s::String) from DynamicQuantities uses eval internally and
# cannot be called during precompilation. Therefore this workload constructs
# ThermoFactory objects WITHOUT string units (units=nothing, output_unit=nothing),
# which falls back to the u"1" macro sentinel and never calls uparse.
# The unit-bearing factory construction happens at runtime via __init__ / THERMO_FACTORIES.

@compile_workload begin
    # Polynomial expression covering the same variables as cp_ft_equation
    # but without string units — avoids uparse during precompilation.
    _f1 = ThermoFactory(:(a₀ + a₁ * T + a₂ / T^2), [:T])
    _sf1 = _f1(; a₀ = 75.0, a₁ = 0.0, a₂ = 0.0)
    _sf1(; T = 298.15)
    _sf1(; T = 500.0)
    _sf1(; T = 298.15, unit = true)

    _f2 = ThermoFactory(:(b₀ * T + b₁ * T^2), [:T])
    _sf2 = _f2(; b₀ = 1.0, b₁ = 0.0)
    _sf2(; T = 298.15)

    # Arithmetic (exercises combine_symbolic and SymbolicFunc binary ops)
    _sf_sum = _sf1 + _sf2
    _sf_sum(; T = 298.15)
    _sf_diff = _sf1 - _sf2
    _sf_diff(; T = 298.15)
end
