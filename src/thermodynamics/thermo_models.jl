"""
    dict_cp_ft_equation

Dictionary of predefined thermodynamic function expressions.

# Contents

  - `:Cp`: Heat capacity expression ``C_p(T)=a_0 + a_1 T + a_2 T^{-2} + a_3 T^{-1/2} + a_4 T^2 + a_5 T^3 + a_6 T^4 + a_7 T^{-3} + a_8 T^{-1} + a_9 T^{1/2} + a_{10}\\ln T``
  - `:CpoverT`: Heat capacity divided by temperature
  - `:H`: Enthalpy expression ``H(T)=\\int C_p(T)\\,\\mathrm{d}T = a_0 T + \\frac{a_1}{2}T^2 - a_2 T^{-1} + 2a_3 T^{1/2} + \\frac{a_4}{3}T^3 + \\frac{a_5}{4}T^4 + \\frac{a_6}{5}T^5 - \\frac{a_7}{2}T^{-2} + a_8\\ln T + \\frac{2}{3}a_9 T^{3/2} + a_{10}T\\ln T - a_{10}T``
  - `:S`: Entropy expression ``S(T)=\\int \\frac{C_p(T)}{T}\\,\\mathrm{d}T = a_0\\ln T + a_1 T - \\frac{a_2}{2}T^{-2} - 2a_3 T^{-1/2} + \\frac{a_4}{2}T^2 + \\frac{a_5}{3}T^3 + \\frac{a_6}{4}T^4 - \\frac{a_7}{3}T^{-3} - a_8 T^{-1} + 2a_9 T^{1/2} + \\frac{a_{10}}{2}(\\ln T)^2``

# Examples

```jldoctest
julia> dict_cp_ft_equation[:Cp]
:(a₀ + a₁ * T + a₂ / T ^ 2 + a₃ / √T + a₄ * T ^ 2 + a₅ * T ^ 3 + a₆ * T ^ 4 + a₇ / T ^ 3 + a₈ / T + a₉ * √T + a₁₀ * log(T))

julia> dict_cp_ft_equation[:CpoverT]
:(a₀ / T + a₁ + a₂ / T ^ 3 + a₃ / T ^ (3 / 2) + a₄ * T + a₅ * T ^ 2 + a₆ * T ^ 3 + a₇ / T ^ 4 + a₈ / T ^ 2 + a₉ / √T + (a₁₀ * log(T)) / T)
```
"""
const dict_cp_ft_equation = Dict(
    :Cp => :(
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
    ),
    :CpoverT => :(
        a₀ / T +
        a₁ +
        a₂ / T^3 +
        a₃ / T^(3 / 2) +
        a₄ * T +
        a₅ * T^2 +
        a₆ * T^3 +
        a₇ / T^4 +
        a₈ / T^2 +
        a₉ / √T +
        a₁₀ * log(T) / T
    ),
    :S => :(
        a₀ * log(T) +
        a₁ * T +
        -a₂ / (2 * T^2) +
        -2 * a₃ / √T +
        (a₄ / 2) * T^2 +
        (a₅ / 3) * T^3 +
        (a₆ / 4) * T^4 +
        -(a₇ /3) / T^3 +
        -a₈ / T +
        2 * a₉ * √T +
        a₁₀ * (log(T))^2 / 2
    ),
    :H => :(
        a₀ * T +
        a₁ * T^2 / 2 +
        - a₂ / T +
        2 * a₃ * √T +
        (a₄ / 3) * T^3 +
        (a₅ / 4) * T^4 +
        (a₆ / 5) * T^5 +
        - (a₇ /2) / T^2 +
        a₈ * log(T) +
        (2 / 3) * a₉ * T^(3 / 2) +
        a₁₀ * T * log(T) - a₁₀ * T
    ),
    :G => :(
        -a₀ * T * log(T) + a₀ * T +
        - (a₁ / 2) * T^2 +
        - (a₂ / 2) / T +
        4 * a₃ * √T +
        - (a₄ / 6) * T^3 +
        - (a₅ / 12) * T^4 +
        - (a₆ / 20) * T^5 +
        - (a₇ / 6) / T^2 +
        a₈ * log(T) +
        - (4 / 3) * a₉ * T^(3 / 2) +
        - a₁₀ * T * (log(T))^2 / 2 + a₁₀ * T * log(T) - a₁₀ * T
    ),
)

"""
    thermo_functions_cp_ft_equation(params, values0; ref=[])

Create a set of standard thermodynamic `ThermoFunction`s from Cp polynomial coefficients.

# Arguments

  - `params`: Coefficient parameters for the Cp polynomial (usually a tuple/array of `a₀..a₁₀`).
  - `values0`: A collection (pairs) of reference values (e.g. `:ΔfH⁰`, `:ΔfG⁰`, `:S⁰`, `:V⁰`).
  - `ref`: Optional reference dictionary specifying reference `:T` and `:P` (for example `[:T => 298.15u"K"]`).

# Returns

  - `Dict` with keys `:Cp⁰`, `:ΔfH⁰`, `:S⁰`, `:ΔfG⁰`, `:V⁰` containing `ThermoFunction` objects.

# Equations

  - Heat capacity (Cp):
      ``C_p(T)=a_0 + a_1 T + a_2 T^{-2} + a_3 T^{-1/2} + a_4 T^2 + a_5 T^3 + a_6 T^4 + a_7 T^{-3} + a_8 T^{-1} + a_9 T^{1/2} + a_{10}\\ln T``
  - Enthalpy (H) obtained by integration:
      ``H(T)=\\int C_p(T)\\,\\mathrm{d}T = a_0 T + \\frac{a_1}{2}T^2 - a_2 T^{-1} + 2a_3 T^{1/2} + \\frac{a_4}{3}T^3 + \\frac{a_5}{4}T^4 + \\frac{a_6}{5}T^5 - \\frac{a_7}{2}T^{-2} + a_8\\ln T + \\frac{2}{3}a_9 T^{3/2} + a_{10}T\\ln T - a_{10}T``
  - Entropy (S) obtained by integrating Cp/T:
      ``S(T)=\\int \\frac{C_p(T)}{T}\\,\\mathrm{d}T = a_0\\ln T + a_1 T - \\frac{a_2}{2}T^{-2} - 2a_3 T^{-1/2} + \\frac{a_4}{2}T^2 + \\frac{a_5}{3}T^3 + \\frac{a_6}{4}T^4 - \\frac{a_7}{3}T^{-3} - a_8 T^{-1} + 2a_9 T^{1/2} + \\frac{a_{10}}{2}(\\ln T)^2``

# The returned formation properties (`ΔfH⁰`, `ΔfG⁰`, `S⁰`) are shifted so they match the provided `values0` at the reference temperature.

Examples

```julia
julia> params = [:a₀ => 210.0u"J/K/mol", :a₁ => 0.12u"J/mol/K^2", :a₂ => -3.07e6u"J*K/mol", :a₃ => 0.0u"J/mol/√K"]
4-element Vector{Pair{Symbol, Quantity{Float64, Dimensions{FRInt32}}}}:
 :a₀ => 210.0 m² kg s⁻² K⁻¹ mol⁻¹
 :a₁ => 0.12 m² kg s⁻² K⁻² mol⁻¹
 :a₂ => -3.07e6 m² kg s⁻² K mol⁻¹
 :a₃ => 0.0 m² kg s⁻² K⁻¹ᐟ² mol⁻¹

julia> values0 = [:Cp⁰ => 210.0u"J/K/mol", :ΔfH⁰ => -2723484.33u"J/mol", :S⁰ => 140u"J/(mol*K)", :ΔfG⁰ => -2480808.197u"J/mol", :V⁰ => 7.84u"J/bar"]
5-element Vector{Pair{Symbol, Quantity{Float64, Dimensions{FRInt32}}}}:
  :Cp⁰ => 210.0 m² kg s⁻² K⁻¹ mol⁻¹
 :ΔfH⁰ => -2.72348433e6 m² kg s⁻² mol⁻¹
   :S⁰ => 140.0 m² kg s⁻² K⁻¹ mol⁻¹
 :ΔfG⁰ => -2.480808197e6 m² kg s⁻² mol⁻¹
   :V⁰ => 7.840000000000001e-5 m³

julia> dtf = thermo_functions_cp_ft_equation(params, values0; ref=[:T => 298.15u"K"])
Dict{Symbol, ThermoFunction{Quantity{Int64, Dimensions{FRInt32}}, F, OrderedCollections.OrderedDict{Symbol, Quantity{Float64, Dimensions{FRInt32}}}} where F} with 5 entries:
  :ΔfH⁰ => -2.80173e6 + 210.0T + 3.07e6 / T + 0.06(T^2) ♢ unit=[m² kg s⁻² mol⁻¹] ♢ ref=[T=298.15 K]
  :S⁰   => -1109.54 + 0.12T + 210.0log(T) + 3.07e6 / (2(T^2)) ♢ unit=[m² kg s⁻² K⁻¹ mol⁻¹] ♢ ref=[T=298.15 K]
  :V⁰   => 7.84e-5 ♢ unit=[m³ mol⁻¹] ♢ ref=[T=298.15 K]
  :ΔfG⁰ => -2.51731e6 + 1319.54T + 1.535e6 / T - 0.06(T^2) - 210.0T*log(T) ♢ unit=[m² kg s⁻² mol⁻¹] ♢ ref=[T=298.15 K]
  :Cp⁰  => 210.0 + 0.12T + -3.07e6 / (T^2) ♢ unit=[m² kg s⁻² K⁻¹ mol⁻¹] ♢ ref=[T=298.15 K]
```
"""
function thermo_functions_cp_ft_equation(params, values0 ; ref=[])
    vars = [:T, :P]
    dict_values0 = Dict(values0)
    dict_ref = Dict(ref)
    Tref = dict_ref[:T]
    T = ThermoFunction(:T; vars=vars, ref=ref)
    with_units = promote_type(typeof.(last.(params))...) <: Quantity

    Cp⁰ = ThermoFunction(dict_cp_ft_equation[:Cp], params; vars=vars, ref=ref)

    H = ThermoFunction(dict_cp_ft_equation[:H], params; vars=vars, ref=ref)
    ΔfH⁰ = H - H(Tref) + dict_values0[:ΔfH⁰]

    S = ThermoFunction(dict_cp_ft_equation[:S], params; vars=vars, ref=ref)
    S⁰ = S - S(Tref) + dict_values0[:S⁰]

    G = ThermoFunction(dict_cp_ft_equation[:G], params; vars=vars, ref=ref)
    ΔfG⁰ = G - G(Tref) + dict_values0[:ΔfG⁰] + (S(Tref) - dict_values0[:S⁰]) * (T - Tref)

    V⁰ = ThermoFunction(:cst => (with_units ? dict_values0[:V⁰] / u"mol" : ustrip(dict_values0[:V⁰])); ref=ref)

    return Dict(:Cp⁰ => Cp⁰, :ΔfH⁰ => ΔfH⁰, :S⁰ => S⁰, :ΔfG⁰ => ΔfG⁰, :V⁰ => V⁰)
end

"""
    thermo_functions_generic_cp_ft(Cpexpr, params, values0; ref=[])

Construct thermodynamic `ThermoFunction`s from a generic Cp expression.

# Arguments

  - `Cpexpr`: An expression (or string) representing the heat capacity functional form in `T`.
  - `params`: Coefficient parameters matching symbols in `Cpexpr`.
  - `values0`: Reference values (pairs) such as `:ΔfH⁰`, `:ΔfG⁰`, `:S⁰`, `:V⁰`.
  - `ref`: Optional reference dictionary specifying `:T` and `:P`.

# Returns

  - `Dict` with keys `:Cp⁰`, `:ΔfH⁰`, `:S⁰`, `:ΔfG⁰`, `:V⁰` where `:Cp⁰` is built from `Cpexpr` and the others are computed by integration
    (with shifts to match `values0` at the reference temperature).

# Equations

  - Given a generic heat capacity expression ``C_p(T)``, the enthalpy and entropy are formed by
      ``H(T)=\\int C_p(T)\\,\\mathrm{d}T\\,,\\qquad S(T)=\\int \\frac{C_p(T)}{T}\\,\\mathrm{d}T\\,,``
      and Gibbs energy
      ``G(T)=-\\int S(T)\\,\\mathrm{d}T\\,,``
      with additive constants chosen so the formation values match `values0` at the reference temperature.

# Examples

```julia
julia> params = [:a₀ => 210.0u"J/K/mol", :a₁ => 0.12u"J/mol/K^2", :a₂ => -3.07e6u"J*K/mol", :a₃ => 0.0u"J/mol/√K"]
4-element Vector{Pair{Symbol, Quantity{Float64, Dimensions{FRInt32}}}}:
 :a₀ => 210.0 m² kg s⁻² K⁻¹ mol⁻¹
 :a₁ => 0.12 m² kg s⁻² K⁻² mol⁻¹
 :a₂ => -3.07e6 m² kg s⁻² K mol⁻¹
 :a₃ => 0.0 m² kg s⁻² K⁻¹ᐟ² mol⁻¹

julia> values0 = [:Cp⁰ => 210.0u"J/K/mol", :ΔfH⁰ => -2723484.33u"J/mol", :S⁰ => 140u"J/(mol*K)", :ΔfG⁰ => -2480808.197u"J/mol", :V⁰ => 7.84u"J/bar"]
5-element Vector{Pair{Symbol, Quantity{Float64, Dimensions{FRInt32}}}}:
  :Cp⁰ => 210.0 m² kg s⁻² K⁻¹ mol⁻¹
 :ΔfH⁰ => -2.72348433e6 m² kg s⁻² mol⁻¹
   :S⁰ => 140.0 m² kg s⁻² K⁻¹ mol⁻¹
 :ΔfG⁰ => -2.480808197e6 m² kg s⁻² mol⁻¹
   :V⁰ => 7.840000000000001e-5 m³

julia> Cpexpr = :(a₀ + a₁ * T + a₂ / T ^ 2 + a₃ / √T + a₄ * T ^ 2 + a₅ * T ^ 3 + a₆ * T ^ 4 + a₇ / T ^ 3 + a₈ / T + a₉ * √T + a₁₀ * log(T))
:(a₀ + a₁ * T + a₂ / T ^ 2 + a₃ / √T + a₄ * T ^ 2 + a₅ * T ^ 3 + a₆ * T ^ 4 + a₇ / T ^ 3 + a₈ / T + a₉ * √T + a₁₀ * log(T))

julia> dtf = thermo_functions_generic_cp_ft(Cpexpr, params, values0; ref=[:T => 298.15u"K"])
Dict{Symbol, ThermoFunction{Quantity{Int64, Dimensions{FRInt32}}, F, OrderedCollections.OrderedDict{Symbol, Quantity{Float64, Dimensions{FRInt32}}}} where F} with 5 entries:
  :ΔfH⁰ => -2.80173e6 + 210.0T + 3.07e6 / T + 0.06(T^2) ♢ unit=[m² kg s⁻² mol⁻¹] ♢ ref=[T=298.15 K]
  :S⁰   => -1109.54 + 0.12T + 210.0log(T) + 1.535e6 / (T^2) ♢ unit=[m² kg s⁻² K⁻¹ mol⁻¹] ♢ ref=[T=298.15 K]
  :V⁰   => 7.84e-5 ♢ unit=[m³ mol⁻¹] ♢ ref=[T=298.15 K]
  :ΔfG⁰ => -2.51731e6 + 1319.54T + 1.535e6 / T - 0.06(T^2) - 210.0T*log(T) ♢ unit=[m² kg s⁻² mol⁻¹] ♢ ref=[T=298.15 K]
  :Cp⁰  => 210.0 + 0.12T + -3.07e6 / (T^2) ♢ unit=[m² kg s⁻² K⁻¹ mol⁻¹] ♢ ref=[T=298.15 K]
```
"""
function thermo_functions_generic_cp_ft(Cpexpr, params, values0 ; ref=[])
    vars = [:T, :P]
    dict_values0 = Dict(values0)
    dict_ref = Dict(ref)
    Tref = dict_ref[:T]
    T = ThermoFunction(:T; vars=vars, ref=ref)
    with_units = promote_type(typeof.(last.(params))...) <: Quantity

    symCpexpr = Num(parse_expr_to_symbolic(Cpexpr, @__MODULE__))
    Cp⁰ = ThermoFunction(symCpexpr, params; vars=vars, ref=ref)

    H = ∫(Cp⁰, :T)
    ΔfH⁰ = H - H(Tref) + dict_values0[:ΔfH⁰]

    CpoverT = ThermoFunction(sum(terms(symCpexpr)/Cp⁰.vars[:T]), params; vars=vars, ref=ref)
    S = ∫(CpoverT, :T)
    S⁰ = S - S(Tref) + dict_values0[:S⁰]

    G = -∫(S, :T)
    ΔfG⁰ = G - G(Tref) + dict_values0[:ΔfG⁰] + (S(Tref) - dict_values0[:S⁰]) * (T - Tref)

    V⁰ = ThermoFunction(:cst => (with_units ? dict_values0[:V⁰] / u"mol" : ustrip(dict_values0[:V⁰])); ref=ref)

    return Dict(:Cp⁰ => Cp⁰, :ΔfH⁰ => ΔfH⁰, :S⁰ => S⁰, :ΔfG⁰ => ΔfG⁰, :V⁰ => V⁰)
end
