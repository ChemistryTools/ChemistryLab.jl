using Revise
using ChemistryLab, Unicode
using DynamicQuantities
using ModelingToolkit
using Plots
using BenchmarkTools

# Construction
expr = :(α + β * T + γ * log(T))
factory = ThermoFactory(expr, [:T])

params = (α=210.0, β=0.0, γ=-3.07e6, T=298.15)

println("=== CONSTRUCTION ===")
@btime $factory(; $params...)  # Premier: ~100 μs, suivants: ~1 ns (cache)

# Exécution
f = factory(; params...)
println("\n=== EXÉCUTION ===")
@btime $f(T=300.0)  # ~5-10 ns (ultra-rapide!)

# Comparaison avec fonction native
f_native(T) = 210.0 + 0.0*T - 3.07e6*log(T)
@btime $f_native(300.0)  # ~5-10 ns (identique!)

# Opérations
f2 = f + 1000
@btime $f2(T=300.0)  # ~5-10 ns (toujours rapide)
