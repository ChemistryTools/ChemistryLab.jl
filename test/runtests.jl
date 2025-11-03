using ChemistryLab
using OrderedCollections
# import Unitful: u, g, cm, K, J, mol, bar, Quantity, uconvert, ustrip, unit, uparse, upreferred, preferunits, @u_str
using DynamicQuantities
using Test

@testset "Construction tests" begin

    @testset "Species tests" begin
        include("species.jl")
    end

    @testset "Formula tests" begin
        include("formulas.jl")
    end

    @testset "Database tests" begin
        include("databases.jl")
    end

    @testset "Parsing tests" begin
        include("parsing_utils.jl")
    end

    # @testset "Reaction tests" begin
    #     include("reactions.jl")
    # end

end