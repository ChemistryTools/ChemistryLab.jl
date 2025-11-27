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
A, atomlist = canonical_stoich_matrix(species; label=:name, display=true) ; # label only for display
A, atomlist = canonical_stoich_matrix(species; label=:symbol, display=true) ;
A, atomlist = canonical_stoich_matrix(species; label=:formula, display=true) ;

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
A, indep_comp = canonical_stoich_matrix(cemspecies; label=:name, display=true) ;
A, indep_comp = canonical_stoich_matrix(cemspecies; label=:symbol, display=true) ;
A, indep_comp = canonical_stoich_matrix(cemspecies; label=:formula, display=true) ;
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
json_file, jls_file = "data/cemdata18-merged.json", "data/cemdata18.jls"
# json_file, jls_file = "data/psinagra-12-07-thermofun.json", "data/psinagra.jls"
try
    global df_substances, df_reactions = deserialize(jls_file)
catch
    # df_elements, df_substances, df_reactions = read_thermofun(json_file; with_units=true, add_species=true, add_reactions=true, all_properties=true, debug=false)
    global df_substances = read_thermofun_substances(json_file; with_units=true, add_species=true, all_properties=true, debug=false)
    global df_reactions = read_thermofun_reactions(json_file, df_substances; with_units=true, add_reactions=true, all_properties=true, debug=false)
    serialize(jls_file, (df_substances, df_reactions))
end
# Construction of Dicts for convenience
dict_species = Dict(zip(df_substances.symbol, df_substances.species))
dict_reactions = Dict(zip(df_reactions.symbol, df_reactions.reaction))
# filter(p->!haskey(p.second, :Cp), dict_species)
# filter(p->haskey(p.second, :Cp) && !iszero(p.second.Cp.a3), dict_species)
df_primaries = extract_primary_species("data/CEMDATA18-31-03-2022-phaseVol.dat")

# Construction of stoich matrix with species from database
given_species = filter(row->row.symbol ∈ split("C3S Portlandite Jennite H2O@"), df_substances)
secondaries = filter(row->row.aggregate_state == "AS_AQUEOUS"
                          && all(k->first(k) ∈ union_atoms(atoms.(given_species.species)), atoms(row.species))
                          && row.symbol ∉ split("H2@ O2@"),
                          df_substances)
all_species = unique(vcat(given_species, secondaries), :symbol)
# species = Species.(all_species.formula)
# candidate_primaries = Species.(df_primaries.formula)
species = [Species(f; symbol=phreeqc_to_unicode(n)) for (f,n) in zip(all_species.formula, all_species.symbol)]
candidate_primaries = [Species(f; symbol=phreeqc_to_unicode(n)) for (f,n) in zip(df_primaries.formula, df_primaries.symbol)]
A, indep_comp, dep_comp = stoich_matrix(species, candidate_primaries; display=true) ;

# Construction of stoich matrix with aqueous species from database
aqueous_species = filter(row->row.aggregate_state == "AS_AQUEOUS", df_substances)
species = [Species(f; symbol=phreeqc_to_unicode(n)) for (f,n) in zip(aqueous_species.formula, aqueous_species.symbol)]
candidate_primaries = [Species(f; symbol=phreeqc_to_unicode(n)) for (f,n) in zip(df_primaries.formula, df_primaries.symbol)]
A, indep_comp, dep_comp = stoich_matrix(species, candidate_primaries; display=true) ;
lr = stoich_matrix_to_reactions(A, indep_comp, dep_comp) ;

# CemSpecies with Sym coef
â, b̂, ĝ = symbols("â b̂ ĝ", real=true)
ox = Dict(:C => â, :S => one(Sym), :A => b̂, :H => ĝ)
CSH = CemSpecies(ox; aggregate_state=AS_CRYSTAL, class=SC_COMPONENT)
numCSH = apply(N, apply(subs, CSH, â=>1.8, b̂=>1, ĝ=>5))
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
â, b̂, ĝ = symbols("â b̂ ĝ", real=true)
CSH = CemSpecies(Dict(:C => â, :S => one(â), :H => ĝ))
C3S = CemSpecies("C3S")
H = CemSpecies("H")
CH = CemSpecies("CH")
r = Reaction([CSH, C3S, H, CH]; equal_sign='→')
 ## application of a function to stoichimetric coefficients (here simplify)
r = apply(simplify, Reaction([C3S, H], [CH, CSH]; equal_sign='→'))
A, _, _ = stoich_matrix([C3S], [CSH, H, CH]; involve_all_atoms=true, display=true) ;
simplify.(A)

# Chen & Brouwers
CSH = CemSpecies(Dict(:C => â, :S => one(â), :A => b̂, :H => ĝ))
HT = CemSpecies("M₅AH₁₃")
HG = CemSpecies("C₆AFS₂H₈")
AFt = CemSpecies("C₆S̄₃H₃₂")
ST = CemSpecies("C₂ASH₈")
AH = CemSpecies("C₄AH₁₃")
A, ox = canonical_stoich_matrix([CSH, HT, HG, AFt, ST, AH]; display=true);
A = typeof(â).(A[1:end-1, 1:end]) # end-1 to remove the line corresponding to water H
oxides = (CemSpecies∘string).(ox[1:end-1])
hydrates = [CSH, HT, HG, AFt, ST, AH]
print_stoich_matrix(A, symbol.(oxides), symbol.(hydrates))
print_stoich_matrix(inv(A), symbol.(hydrates), symbol.(oxides))
Mhyd = ustrip.(getproperty.(hydrates, :M))
Mox = ustrip.(getproperty.(oxides, :M))
B = Mox .* A .* inv.(Mhyd)'
# or directly
B, ox = canonical_stoich_matrix([CSH, HT, HG, AFt, ST, AH]; mass=true, display=true) ;
B = B[1:end-1, 1:end] # to remove the H line
print_stoich_matrix(B, "m_" .* symbol.(oxides), "m_" .* symbol.(hydrates))
print_stoich_matrix(subs.(inv(B), Ref(Dict(â=>1.8, b̂=>1, ĝ=>4))), "m_" .* symbol.(hydrates), "m_" .* symbol.(oxides))

# Alkane combustion with SymPy
n = symbols("n", real=true) ;
CₙH₂ₙ₊₂ = Species(:C => n, :H => 2n+2) ;
O₂ = Species("O₂") ;
H₂O = Species("H₂O") ;
CO₂ = Species("CO₂") ;
r = Reaction([CₙH₂ₙ₊₂, O₂], [H₂O, CO₂])
apply(factor, r)
println(2r)
for vn in 1:9 println("n=$vn ⇒ ", apply(subs, r, n=>vn)) end
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
for vn in 1:9 println("n=$vn ⇒ ", apply(substitute, r, n=>vn)) end
for vn in 1:9 println("n=$vn ⇒ ", apply(stoich_coef_round∘(x->x isa Num ? x.val : x)∘substitute, r, n=>vn)) end

# Example from https://github.com/thermohub/chemicalfun
formulas = ["Ca+2", "Fe+2", "Fe|3|+3", "H+", "OH-", "SO4-2", "CaSO4@", "CaOH+", "FeO@", "HFe|3|O2@", "FeOH+", "Fe|3|OH+2", "H2O@",  "FeS|-2|", "FeS|0|S|-2|", "S|4|O2"] ;
species = Species.(formulas) ;
candidate_primaries = species[1:6] ;
A, indep_comp, dep_comp = stoich_matrix(species; display=true) ;
B, indep_comp, dep_comp = stoich_matrix(species; mass=true, display=true) ;
lr = stoich_matrix_to_reactions(A, indep_comp, dep_comp) ;

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
            println(s, "  ∙  ", row.symbol)
        end
    end
end

for row in eachrow(df_reactions)
    re = row.reaction
    for s in keys(re.products)
        if !haskey(s, :Cp)
            println(s, "  ∙  ", row.symbol)
        end
    end
end

# Check consistency of logKr at Tref and logKr_Tref of the database
for r in df_reactions.reaction
    println(r, " → ", r.logKr(), " == ", r.logKr_Tref)
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
        println(collect(keys(r))[1].symbol, " ", r)
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
A, indep_comp = canonical_stoich_matrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻]; display=true) ;
A, indep_comp, dep_comp = stoich_matrix([H₂O, H⁺, OH⁻, CO₂, HCO₃⁻, CO₃²⁻]; display=true) ;
