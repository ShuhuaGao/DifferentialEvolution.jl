using DifferentialEvolution
using Test

@testset "DifferentialEvolution.jl" begin
    include("nsga_ii_test.jl")
    include("constraint_test.jl")
end
