# Advanced Topics

This page covers advanced usage patterns and techniques for working with ChemistryLab: complex formula transformations, cement species, reaction algebra, and programmatic database operations.

## Advanced Formula Manipulation

### Working with Quantities (units)

If your project uses `DynamicQuantities` or `Unitful`, you can attach units to stoichiometric coefficients and apply transformations that preserve dimensions:

```julia
using ChemistryLab, DynamicQuantities

# Create formulas with quantities
comp_with_units = OrderedDict(:H => 2.0u"g/mol", :O => 1.0u"g/mol")
f = Formula(comp_with_units)

# Apply a function across coefficients (units preserved where possible)
f_scaled = apply(x -> x * 2, f)
```

### Arithmetic with fractional stoichiometry

ChemistryLab preserves rational coefficients when parsing fractional formulas and when doing arithmetic:

```julia
f = Formula("H1/2O")
composition(f)  # OrderedDict(:H => 1//2, :O => 1)

f2 = f * 2
composition(f2) # OrderedDict(:H => 1, :O => 2)  — now simplified to integers
```

### Converting between notations

Switch between Phreeqc (plain text), Unicode (pretty), and colored output:

```julia
f = Formula("Fe+3")

expr(f)       # "Fe+3"
phreeqc(f)    # "Fe+3"
unicode(f)    # "Fe³⁺"
colored(f)    # colored terminal output (see it in REPL)

# Convert between notations programmatically
phreeqc_to_unicode("SO4-2")  # "SO₄²⁻"
unicode_to_phreeqc("SO₄²⁻")  # "SO4-2"
```

## Advanced Species Operations

### CemSpecies and oxide-to-atom decomposition

`CemSpecies` represents a species in cement nomenclature (as oxides, e.g., C3S = 3CaO·SiO2) and automatically converts to/from atomic composition:

```julia
# Create a cement species from oxide formula
c3s = CemSpecies("C3S")
oxides(c3s)       # OrderedDict(:C => 3, :S => 1)  — cement components
atoms(c3s)        # OrderedDict(:Ca => 3, :Si => 1, :O => 8)  — atomic decomposition

# Convert between Species and CemSpecies
s = Species("Ca3SiO8")
cem_s = CemSpecies(s)  # automatically decomposes to oxides
```

### Type conversion and promotion

Convert numeric types in a Species and use promotion to work with mixed species types:

```julia
s1 = Species{Int}("H2O")
s2 = Species{Float64}(s1)  # convert coefficients to Float64

# Promotion rules allow mixed arithmetic
s_aqueous = Species("H2O", aggregate_state=AS_AQUEOUS)
s_crystal = Species("H2O", aggregate_state=AS_CRYSTAL)
# Both are Species and can be used together
```

### Custom properties and thermodynamic data

Attach arbitrary properties to species (molar mass is auto-calculated; add more):

```julia
water = Species("H2O")
water[:Cp] = 75.3  # J/(mol·K)
water[:ΔfH0] = -285.8  # kJ/mol
water[:S0] = 69.9   # J/(mol·K)

# Access via bracket or dot notation
water[:Cp]   # 75.3
water.Cp     # 75.3 (same)

# Check if a property exists
haskey(water, :Cp)  # true
```

## Advanced Reaction Operations

### Reaction algebra (arithmetic)

Combine and transform reactions using operator overloading:

```julia
r1 = Reaction("H2O = H+ + OH-")
r2 = Reaction("H2O = H2 + 1/2 O2")

r_sum = r1 + r2      # combine all species and coefficients
r_diff = r1 - r2     # subtract r2 from r1 (reverse coefficients of r2)

# Scalar multiplication
r_scaled = 2 * r1    # double all coefficients
r_half = r1 / 2      # halve all coefficients
```

### Reaction simplification

Eliminate species appearing on both sides (cancel them):

```julia
# Build a reaction with redundant species
reac = OrderedDict(Species("A") => -1, Species("B") => -1)
prod = OrderedDict(Species("B") => 1, Species("C") => 1)
r = Reaction(reac, prod)
# r now has: A + B = B + C

# Simplify
r_simple = simplify_reaction(r)
# r_simple now has: A = C  (B cancelled)
```

### Building reactions from species lists

Construct a reaction by specifying independent and dependent species (the algorithm will solve for stoichiometry):

```julia
# Water splitting reaction
h2o = Species("H2O")
h2 = Species("H2")
o2 = Species("O2")

# Declare reactants and products; stoichiometry is solved automatically
r = Reaction([h2o], [h2, o2]; auto_scale=true)
# Result: 2H2O = 2H2 + O2  (stoichiometry auto-scaled to integers)
```

## Advanced Stoichiometric Matrix Operations

### Mass-based stoichiometric matrices

By default, stoichiometric matrices use atom counts. Set `mass=true` to get mass-balanced matrices (coefficients adjusted by molar masses):

```julia
species = [Species("H2O"), Species("H2"), Species("O2")]
A_atoms, _, _ = stoich_matrix(species; display=false, mass=false)
A_mass, _, _ = stoich_matrix(species; display=false, mass=true)
# A_mass rows scaled by atomic masses — different physical meaning
```

### Extracting independent and dependent species

The stoichiometric matrix algorithm automatically identifies independent (basis) and dependent (derived) species:

```julia
species = [Species("Ca2+"), Species("OH-"), Species("Ca(OH)2")]
A, indep, dep = stoich_matrix(species; display=false)

# indep: basis species (linearly independent)
# dep: species expressed as combinations of indep
# A: matrix where A[i,j] = coefficient of indep[i] in dep[j]
```

### Reordering and redox handling

The algorithm handles charged species and redox reactions; use `involve_all_atoms=true` to include all atoms in the matrix (default restricts to atoms in the primary species):

```julia
# Redox reaction with electron transfer
species = [Species("Fe2+"), Species("Fe3+"), Species("e-")]
A, indep, dep = stoich_matrix(species; display=false, involve_all_atoms=true)
# Includes :Zz (charge marker) and :e- (electrons)
```

## Advanced Database Operations

### Reading and filtering ThermoFun data

Load a ThermoFun JSON file and filter substances or reactions by criteria:

```julia
# Read ThermoFun JSON
df_elements, df_substances, df_reactions = read_thermofun("data/cemdata18-thermofun.json"; with_units=false, debug=false)

# Filter aqueous substances
aqueous = filter(row -> row.aggregate_state == "AS_AQUEOUS", df_substances)

# Filter by charge
charged_species = filter(row -> row.charge != 0, df_substances)
```

### Merging custom reactions into a ThermoFun source

Build a stoichiometric matrix, convert to reactions, and merge into an existing ThermoFun JSON:

```julia
# Define new species
new_species = [Species("CustomA"), Species("CustomB"), Species("CustomC")]
A, indep, dep = stoich_matrix(new_species; display=false)

# Convert to reactions and merge (pseudo-code; exact merge API varies)
# new_reactions = stoich_matrix_to_reactions(A, indep, dep; display=false)
# merged_json = merge_reactions(original_json, new_reactions)
```

### Cemdata .dat parsing and extraction

Extract primary species and phases from a Phreeqc / Cemdata .dat file:

```julia
# Extract primary aqueous species
df_primary = extract_primary_species("path/to/file.dat")

# Parse phase equilibria (PHASES section)
phases = parse_phases(dat_content)
```

## Tips for complex workflows

1. **Batch operations on Species/Formulas**: Use `apply(func, formula)` or `apply(func, species)` to transform stoichiometric values or properties in bulk.

2. **Debugging equations**: Use `display=true` (default) in `stoich_matrix` and related functions to print colored, formatted matrices and equations to the REPL.

3. **Type stability**: Prefer homogeneous numeric types (`Species{Float64}` or `Species{Rational}`) across collections for performance.

4. **Combining databases**: Use `read_thermofun` to load multiple sources, then concatenate DataFrames or merge reactions as needed.

5. **Custom properties as metadata**: Attach source information, uncertainty, or application-specific data as custom properties — they're preserved during conversions and copies.

## Next steps

- Explore the reference documentation for full API signatures and options.
- See `docs/src/example/` for complete worked examples (cement hydration calculations, equilibrium problems, etc.).
- If you encounter edge cases or need custom workflows, open an issue or extend the examples here.
