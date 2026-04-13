using ModelingToolkit
using SymPy

Cpexpr = :(
    a₀ +
        a₁ * T +
        a₂ / T^2 +
        a₃ / √T +
        a₄ * T^2 +
        a₅ * T^3 +
        a₆ * T^4 +
        a₇ / T^3 +
        a₈ / T +
        a₉ * √T +
        a₁₀ * log(T)
)

Cp = Num(parse_expr_to_symbolic(Cpexpr, @__MODULE__))

T = Num(parse_expr_to_symbolic(:T, @__MODULE__))

int = expr -> symbolics_to_sympy(sympy_integrate(expr, T))

H = int(Cp)

S = int(sympy_simplify(Cp / T))

G = -int(S)
