using Revise, ChemistryLab, Unicode
using DynamicQuantities
using ModelingToolkit
using Serialization
import SymPy: symbols, Sym, N, subs, factor, simplify

# Formula
fgen = Formula("A1//2B3C0.4")
convert(Float64, fgen)
ATOMIC_ORDER # provides the order of atoms in formulas

# Species
H₂O = Species("H₂O"; name="Water", symbol="H₂O@", aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)
HSO₄⁻ = Species("HSO₄⁻"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE)
CO₂ = Species(Dict(:C=>1, :O=>2); name="Carbon dioxide", symbol="CO₂⤴", aggregate_state=AS_GAS, class=SC_GASFLUID)
species = [H₂O, HSO₄⁻, CO₂] ;
CSM = CanonicalStoichMatrix(species)
pprint(CSM; label=:name)
pprint(CSM; label=:symbol)
pprint(CSM; label=:formula)

water_without_name_symbol = Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT)
water_without_name_symbol == H₂O # true since atoms, aggregate_state and class are equal despite instances are different
vapour = Species("H₂O"; name="Vapour", symbol="H₂O⤴", aggregate_state=AS_GAS, class=SC_GASFLUID)
vapour == H₂O # false since aggregate_state or class are different despite atoms are identical

# CemSpecies
OXIDE_ORDER # provides the order of oxides in cement formulas
C3S = CemSpecies("C3S"; name="Alite", symbol="C₃S", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
C2S = CemSpecies("C₂S"; name="Belite", symbol="C₂S", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
C3A = CemSpecies("C3A"; name="Aluminate", symbol="C₃A", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
C4AF = CemSpecies(Dict(:C=>4, :A=>1, :F=>1); name="Ferrite", symbol="C₄AF", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
cemspecies = [C3S, C2S, C3A, C4AF]
CSM = CanonicalStoichMatrix(cemspecies)
pprint(CSM; label=:name)
pprint(CSM; label=:symbol)
pprint(CSM; label=:formula)
 # conversion CemSpecies → Species always possible
spC3S = Species(C3S)
spC3S == C3S # true since atoms, aggregate_state and class are identical
 # conversion Species → CemSpecies possible only if the species can be decomposed in cement oxides
spCH = Species("Ca(OH)2"; name="Portlandite", symbol="CH", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
CH = CemSpecies(spCH)
spCH == CH # true
spCH == CemSpecies("CH") # false since aggregate_state and class are undef
spCH == CemSpecies("CH"; aggregate_state=AS_CRYSTAL, class=SC_COMPONENT) # true even though names and/or symbols do not coincide
try CemSpecies(Species("Ca(OH)")) catch; "ERROR: Ca(OH) cannot be decomposed in cement oxides" end
CemSpecies(Species("CaCO3"; name="Calcite", aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)) # ok here

# Thermofun input
function extract_database(json_file)
    jls_file = splitext(json_file)[1] * ".jls"
    df_substances, df_reactions = nothing, nothing
    try
        df_substances, df_reactions = deserialize(jls_file)
    catch
        # df_elements, df_substances, df_reactions = read_thermofun(json_file; with_units=true, add_species=true, add_reactions=true, all_properties=true, debug=false)
        df_substances = read_thermofun_substances(json_file; with_units=true, add_species=true, all_properties=true, debug=false)
        df_reactions = read_thermofun_reactions(json_file, df_substances; with_units=true, add_reactions=true, all_properties=true, debug=false)
        serialize(jls_file, (df_substances, df_reactions))
    end
    dict_species = Dict(zip(df_substances.symbol, df_substances.species))
    dict_reactions = Dict(zip(df_reactions.symbol, df_reactions.reaction))
    return df_substances, df_reactions, dict_species, dict_reactions
end
df_substances, df_reactions, dict_species, dict_reactions = extract_database("data/cemdata18-merged.json")
df_substances_psi, df_reactions_psi, dict_species_psi, dict_reactions_psi = extract_database("data/psinagra-12-07-thermofun.json")
df_substances_aq17, df_reactions_aq17, dict_species_aq17, dict_reactions_aq17 = extract_database("data/aq17-thermofun.json")

# Extraction of primaries from .dat
df_primaries = extract_primary_species("data/CEMDATA18-31-03-2022-phaseVol.dat")

# Construction of stoich matrix with species from database
given_species = filter(row->row.symbol ∈ split("C3S Portlandite Jennite H2O@"), df_substances)
secondaries = filter(row->row.aggregate_state == "AS_AQUEOUS"
                          && all(k->first(k) ∈ union_atoms(atoms.(given_species.species)), atoms(row.species))
                          && row.symbol ∉ split("H2@ O2@"),
                          df_substances)
species = unique(vcat(given_species, secondaries), :symbol).species
candidate_primaries = [s == "Zz" ? Species("Zz") : dict_species[s] for s in df_primaries.symbol]
SM = StoichMatrix(species, candidate_primaries)
pprint(SM)
list_reactions = reactions(SM)
pprint(list_reactions)

# Construction of stoich matrix with aqueous species from database
aqueous_species = filter(row->row.aggregate_state == "AS_AQUEOUS", df_substances)
species = aqueous_species.species
candidate_primaries = [s == "Zz" ? Species("Zz") : dict_species[s] for s in df_primaries.symbol]
SM = StoichMatrix(species, candidate_primaries)
pprint(SM)
list_reactions = reactions(SM)
pprint(list_reactions)

# CemSpecies with Sym coef
a, b, g = symbols("a b g", real=true)
ox = Dict(:C => a, :S => one(Sym), :A => b, :H => g)
CSH = CemSpecies(ox; aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
numCSH = apply(N, apply(subs, CSH, a=>1.8, b=>1, g=>5))
floatCSH = apply(x->convert(Float64, x), numCSH) # only coefficients of oxides are converted to Float64 here not those of atoms

# Conversion to cement notation with species from database
CemSpecies(H₂O)
CemSpecies(vapour)
df_Jennite = filter(row->row.symbol == "Jennite", df_substances)
Jennite = Species(df_Jennite.formula[1]; name=df_Jennite.name[1], symbol=df_Jennite.symbol[1], aggregate_state=eval(Meta.parse(df_Jennite.aggregate_state[1])), class=eval(Meta.parse(df_Jennite.class[1])))
cemJennite = CemSpecies(Jennite)
Jennite == cemJennite

# Extraction of properties
cemJennite.ΔfG⁰ = df_Jennite.ΔfG⁰[1].values*uparse(df_Jennite.ΔfG⁰[1].units)
cemJennite

# Equation parsing
equation = "13H⁺ + NO₃⁻ + CO₃²⁻ + 10e⁻ = 6H₂O@ + HCN@"
 ## simple parsing giving Dicts
reac, prod, equal_sign = parse_equation(equation)
 ## construction of a Reaction struct from a string
r = Reaction(equation)
 ## simplification of a Reaction
r = Reaction("2H₂O → H⁺ + OH⁻ + H₂O")
println(simplify_reaction(r))
 ## construction of a Reaction of only CemSpecies by CemReaction
eqC3S = "C₃S + 5.3H = 1.3CH + C₁.₇SH₄"
rC3S = CemReaction(eqC3S)
 ## construction of a Reaction by operations on Species
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
CSH = CemSpecies("C1.7SH4")
r = C3S + 5.3H ↔ 1.3CH + CSH
 ## construction of a Reaction by a balance calculation
r = Reaction([C3S, H, CH, CSH]; equal_sign='→')
Reaction(CemSpecies.(["C3S", "H", "CH", "C1.8SH4"]))
for c_over_s in 1.5:0.1:2.
    println(Reaction(CemSpecies.(["C3S", "H", "CH", "C$(c_over_s)SH4"])))
end
 ## construction of a Reaction by a balance calculation with symbolic numbers
a, b, g = symbols("a b g", real=true)
CSH = CemSpecies(Dict(:C => a, :S => one(a), :H => g))
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
r = Reaction([CSH, C3S, H, CH]; equal_sign='→')
 ## application of a function to stoichimetric coefficients (here simplify)
r = apply(simplify, Reaction([C3S, H], [CH, CSH]; equal_sign='→'))
SM = StoichMatrix([C3S], [CSH, H, CH]; involve_all_atoms=true) ; pprint(SM)
pprint(apply(x->simplify.(x), SM))

# Chen & Brouwers
CSH = CemSpecies(Dict(:C => a, :S => one(a), :A => b, :H => g))
HT = CemSpecies("M₅AH₁₃")
HG = CemSpecies("C₆AFS₂H₈")
AFt = CemSpecies("C₆S̄₃H₃₂")
ST = CemSpecies("C₂ASH₈")
AH = CemSpecies("C₄AH₁₃")
CSM = CanonicalStoichMatrix([CSH, HT, HG, AFt, ST, AH]); pprint(CSM)
A, ox = CSM.A, CSM.primaries
A = typeof(a).(A[1:end-1, 1:end]) # end-1 to remove the line corresponding to water H
oxides = (CemSpecies∘string).(ox[1:end-1])
hydrates = [CSH, HT, HG, AFt, ST, AH]
pprint(A, symbol.(oxides), symbol.(hydrates))
pprint(inv(A), symbol.(hydrates), symbol.(oxides))
Mhyd = ustrip.(us"g/mol", getproperty.(hydrates, :M))
Mox = ustrip.(us"g/mol", getproperty.(oxides, :M))
B = Mox .* A .* inv.(Mhyd)'
# or directly
mCSM = mass_matrix(CSM)
B, ox = mCSM.A, mCSM.primaries
B = B[1:end-1, 1:end] # to remove the H line
pprint(B, "m_" .* symbol.(oxides), "m_" .* symbol.(hydrates))
pprint(subs.(inv(B), Ref(Dict(a=>1.8, b=>1, g=>4))), "m_" .* symbol.(hydrates), "m_" .* symbol.(oxides))

# Alkane combustion with SymPy
n = symbols("n", real=true) ;
CₙH₂ₙ₊₂ = Species(:C => n, :H => 2n+2) ;
O₂ = Species("O₂") ;
H₂O = Species("H₂O") ;
CO₂ = Species("CO₂") ;
r = Reaction([CₙH₂ₙ₊₂, O₂], [H₂O, CO₂])
apply(factor, r)
pprint(2r)
for vn in 1:9 print("n=$vn ⇒ "); println(colored(apply(subs, r, n=>vn))) end
println(Reaction([CₙH₂ₙ₊₂, O₂], [H₂O, CO₂]; side=:products))
println(Reaction([CₙH₂ₙ₊₂, O₂], [H₂O, CO₂]; side=:reactants))
@show r[O₂]

# Alkane combustion with Symbolics
@variables n
CₙH₂ₙ₊₂ = Species(:C => n, :H => 2n+2)
O₂ = Species("O₂")
H₂O = Species("H₂O")
CO₂ = Species("CO₂")
r = Reaction([CₙH₂ₙ₊₂, O₂], [H₂O, CO₂])
for vn in 1:9 print("n=$vn ⇒ "); println(colored(apply(substitute, r, n=>vn))) end
for vn in 1:9 print("n=$vn ⇒ "); println(colored(apply(stoich_coef_round∘(x->x isa Num ? x.val : x)∘substitute, r, n=>vn))) end

# Example from https://github.com/thermohub/chemicalfun
formulas = ["Ca+2", "Fe+2", "Fe|3|+3", "H+", "OH-", "SO4-2", "CaSO4@", "CaOH+", "FeO@", "HFe|3|O2@", "FeOH+", "Fe|3|OH+2", "H2O@",  "FeS|-2|", "FeS|0|S|-2|", "S|4|O2"] ;
species = Species.(formulas) ;
candidate_primaries = species[1:6] ;
SM = StoichMatrix(species)
pprint(SM)
list_reactions = reactions(SM)
pprint(list_reactions)


# Callable
 # with units (coefficient units should be consistent with the basis of functions provided in thermofun database)
coeffs = [:a₀ => 210.0u"J/K/mol", :a₁ => 0.12u"J/mol/K^2", :a₂ => -3.07e6u"J*K/mol", :a₃ => 0.0u"J/mol/√K"]
cemJennite.Cp = ThermoFunction(:Cp, coeffs; ref=[:T=>298.15u"K", :P=>1u"bar"])
@show cemJennite.Cp ;
@show cemJennite.Cp(298.15u"K") ;
@show cemJennite.Cp() ; # application by default on Tref
 # same without units
cemJennite.Cp = ThermoFunction(:Cp, ustrip.(coeffs); ref=[:T=>298.15u"K", :P=>1u"bar"])
@show cemJennite.Cp ;
@show cemJennite.Cp(298.15) ;
@show cemJennite.Cp() ; # application by default on Tref

using Plots
lT = ((0:1:100) .+ 273.15).*u"K"
@time plot(ustrip.(lT), ustrip.(dict_species["Jennite"].ΔfH⁰.(lT)))

# Check which species involved in reactions have not been previously constructed in the list of substances (in this case they are built on-the-fly and don't have thermo properties)
for row in eachrow(df_reactions)
    re = row.reaction
    for s in keys(re.reactants)
        if !haskey(s, :Cp)
            println(colored(s), "  ∙  ", row.symbol)
        end
    end
end

for row in eachrow(df_reactions)
    re = row.reaction
    for s in keys(re.products)
        if !haskey(s, :Cp)
            println(colored(s), "  ∙  ", row.symbol)
        end
    end
end

# Check consistency of logKr at Tref and logKr_Tref of the database
for r in df_reactions.reaction
    println(colored(r), " → ", r.logKr(), " == ", r.logKr_Tref)
end

coeffs = [:a₀ => 210.0u"J/K/mol", :a₁ => 0.12u"J/mol/K^2", :a₂ => -3.07e6u"J*K/mol", :a₃ => 0.0u"J/mol/√K"]
Cp = ThermoFunction(:Cp, coeffs)

rate = ThermoFunction(:((c₁+c₂*t)/(c₃+c₄*√t)), [:c₁ => 1.0, :c₂ => 2.0u"1/s", :c₃ => 3.0, :c₄ => 4.0u"1/√s"])
rate(1u"s")

r = dict_reactions["Cal"]
Tref = 298.15u"K"
ΔrG⁰(T) = sum(ν*s.ΔfG⁰(T) for (s,ν) in r)
logK = -ΔrG⁰(Tref)/(Constants.R*Tref)/log(10)
K(T) = 10^(-ΔrG⁰(T)/(Constants.R*T)/log(10))
pK(T) = ΔrG⁰(T)/(Constants.R*T)/log(10)
r.logKr()
plot(ustrip.(lT), ustrip.(ΔrG⁰.(lT)))
for (re,ν) in r.reactants
    println(re, " ΔfG⁰=", re.ΔfG⁰(Tref))
end
for (pr,ν) in r.products
    println(pr, " ΔfG⁰=", pr.ΔfG⁰(Tref))
end
for (s, ν) in r
    println(s, " ΔfG⁰=", s.ΔfG⁰(Tref))
end

for r in df_reactions.reaction
    if true # abs((r.logKr_Tref -r.logKr())/r.logKr_Tref)>0.01
        println(collect(keys(r))[1].symbol)
        println(colored(r))
        println("     → logKr given at $(Tref) = ", r.logKr_Tref)
        println("     → logKr calculated at $(Tref) = ", r.logKr())
        try
            println("     → logKr=-ΔrG⁰(Tref)/(R Tref ln(10)) = ", -r.ΔrG⁰(Tref)/(Constants.R*Tref)/log(10))
        catch
            println("     → logKr=-ΔrG⁰(Tref)/(R Tref ln(10)) = XXXXXXXXXXXXXXXXXXXXX")
        end
    end
end

# Construction of a Reaction from a string and a list of species (possibly with thermodynamic properties from a database)
r = Reaction("NaOH → Na+ + OH-"; species_list = df_substances.species)
for re in keys(r.reactants)
    display(re)
end
for pr in keys(r.products)
    display(pr)
end

# Example from Leal et al. 2017
H₂O = dict_species["H2O@"]
H⁺ = dict_species["H+"]
OH⁻ = dict_species["OH-"]
CO₂ = dict_species["CO2@"]
HCO₃⁻ = dict_species["HCO3-"]
CO₃²⁻ = dict_species["CO3-2"]
for sp in (:H₂O, :H⁺, :OH⁻, :CO₂, :HCO₃⁻, :CO₃²⁻)
    symsp = String(sp)
    @eval $sp = Species($(String(sp)))
end
CSM = CanonicalStoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻]); pprint(CSM)
SM = StoichMatrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻]); pprint(SM)

# Example from Miron et al. 2023, https://doi.org/10.21105/joss.04624
p1 = plot(xlabel="Temperature [°C]", ylabel="Cp⁰ [J K⁻¹]", xticks=0:50:250, yticks=-2000:200:2000)
for sp ∈ getindex.(Ref(dict_species_aq17), split("Na+ Ca+2 SiO2@ CO3-2 OH-"))
    plot!(p1, θ->sp.Cp⁰(273.15+θ), 0:0.1:250, label=unicode(sp))
end
p2 = plot(xlabel="Temperature [°C]", ylabel="log₁₀K⁰", xticks=0:50:250, yticks=-20:2:10)
for lsp ∈ [
            "Calcite Ca+2 CO3-2",
            "H2O@ H+ OH-",
            "NaCl@ Na+ Cl-",
            "Al+3 H2O@ AlOH+2 H+",
          ]
    rr = Reaction(getindex.(Ref(dict_species_aq17), split(lsp)))
    plot!(p2, θ->rr.logK⁰(273.15+θ), 0:0.1:250, label=rr.equation)
end
plot(p1, p2, layout = (1, 2))
