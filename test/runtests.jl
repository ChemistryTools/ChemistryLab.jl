using ChemistryLab
using DynamicQuantities
using SymbolicNumericIntegration
using Test
using TimerOutputs

macro testsection(str, block)
    return quote
        @timeit "$($(esc(str)))" begin
            @testset "$($(esc(str)))" begin
                $(esc(block))
            end
        end
    end
end

reset_timer!()

@testsection "Construction tests" begin

    include("species.jl")
    include("formulas.jl")
    include("stoich_matrices.jl")
    include("databases.jl")
    include("parsing_utils.jl")
    include("reactions.jl")

end

@testsection "Speciation tests" begin
    include("speciation.jl")
end

@testsection "System and State tests" begin
    include("chemical_systems.jl")
    include("chemical_states.jl")
end

@testsection "Thermodynamics tests" begin
    include("thermodynamics.jl")
    include("hkf.jl")
end

@testsection "Equilibrium tests" begin
    include("activities.jl")
    include("solid_solutions.jl")
end

@testsection "Utils tests" begin
    include("utils.jl")
end

print_timer()
println()
