module DifferentialEvolution
using ProgressMeter

# Write your package code here.
export DEEvolver, pop, fitness, senses, nindividuals, nobjectives, best_fitness, best_individual,
    StdDE, StdOptions, mutate_rand_1, mutate_rand_1!, crossover_binomial, StdSoController, evolve!,
    nondominated_sort, get_nondomination_rank, crowded_compare, dominate, assign_crowding_distance,
    bounce_back!, demo!

include("common.jl")
include("utils.jl")
include("std_de.jl")
include("mutation.jl")
include("crossover.jl")
include("constraint.jl")
include("mo/nsga_ii.jl")
include("mo/mo.jl")
include("mo/demo.jl")
include("mo/debi.jl")
end
