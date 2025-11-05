using ChemistryLab
using DynamicQuantities
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
    include("stoich_matrix.jl")
    include("databases.jl")
    include("parsing_utils.jl")
    include("reactions.jl")

end

print_timer()
println()



