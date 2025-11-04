using Symbolics
using DynamicQuantities

const g, cm, K, J, mol, bar = us"g", us"cm", us"K", us"J", us"mol", us"bar"
DynamicQuantities.uconvert(u, ::Missing) = missing

struct ThermoFunction{N,F}
    symbolic::N
    func::F
end

degrees = [0, 1, -2, -0.5, 2, 3, 4, -3, -1, 0.5, :log]
coeffs = [210.0J/K/mol, 0.120J/mol/K^2, -3.07e6J*K/mol, 0.0J/mol/√K]

function ThermoFunction()
    @variables T, a₀, a₁, a₂
    symbolic = a₀ + a₁ * T + a₂ / T^2
    func = eval(build_function(symbolic, T, a₀, a₁, a₂; expression=Val(true)))
    return ThermoFunction(symbolic, func)
end

function deriv(expr, var)
    return Symbolics.derivative(expr, (@eval @variables var)[1])
end

function derivT(expr)
    @variables T
    return Symbolics.derivative(expr, T)
end

function (tf::ThermoFunction)(args...)
    return tf.func(args...)
end
