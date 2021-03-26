#===================================================================================================
Standard differential evolution
===================================================================================================#
import Base
using Random: shuffle!
using Formatting
using Statistics: std, mean
using StaticArrays: MVector

"""
    StdOptions{T<:AF}

Standard options with two parameters, `F` (scaling) and `Cr` (crossover rate).
This type can be constructed with either positional arguments or keyword arguments.
"""
Base.@kwdef mutable struct StdOptions{T<:AF} <: DEOptions
    F::T
    Cr::T
end


"""
    StdDE{TP<:AM{<:AF}, TF<:AVoM{<:Real}, TS<:ASoV{<:Integer}, TO<:DEOptions} <: DEEvolver

The standard DE type with the following fields:
- `pop::TP`: a matrix-type population, each column of which is an individual (vector) in DE
- `fitness::TF`: fitness of each individual. For single-objective DE, `fitness` is a vector-like 
    structure; for multi-objective DE, `fitness` is of a matrix type, whose column contains values 
    of multiple objectives for each individual.
- `senses::TS`: integer for a single objective or vector of integers for multiple objectives. 
    An individual's strength is defined by `fitness` times `senses`. Only values 1 (maximization) 
    and -1 (minimization) are allowed. To perform minimization of an 
    objective, set a negative sense, i.e., -1.
- `bounds::TP`: a `d-by-2` matrix for a `d`-dim DE, each row containing the lower and upper boundaries
    of a variable.
- `options::TO`: options which usually contain the control parameters. See also [`StdOptions`](@ref).

In practice, more convenient types are usually used like [`StdSoDE`](@ref) for single-objective DE 
and [`StdMoDE`](@ref) for multiobjective DE.
"""
struct StdDE{TP<:AM{<:AF}, TF<:AVoM{<:Real}, TS<:ASoV{<:Integer}, 
        TO<:DEOptions} <: DEEvolver
    pop::TP
    fitness::TF
    senses::TS
    bounds::TP
    options::TO
    _rand_indices::Vector{Int}  # to avoid frequent memory allocation
    function StdDE(pop::AM{TV}, fitness::AVoM{<:Real}, senses::ASoV{<:Integer}, bounds::AM{TV}, 
            options::DEOptions) where {TV<:AF}
        n = size(pop, 2)
        @assert n >= 5  "At least 5 individuals"
        @assert ndims(fitness) <= 2
        @assert n == (ndims(fitness) == 1 ? length(fitness) : size(fitness, 2))
        @assert size(pop, 1) == size(bounds, 1)
        @assert all(s -> s == 1 || s == -1, senses) "A sense can only be 1 or -1"
        if senses isa Integer
            @assert fitness isa AbstractVector
        else 
            @assert fitness isa AbstractMatrix
        end
        new{typeof(pop), typeof(fitness), typeof(senses), typeof(options)}(pop, fitness, senses, 
            bounds, options, collect(1:n))
    end
end

"""
    StdDE(::Type{TF}, nindividuals::Int, bounds::AM{<:AF}, opt::DEOptions; weights::ASoV{<:Real}) where {TF<:Real}

Create a standard DE evolver with a population of `nindividuals` individuals and options specified by 
`opt`. The length of each individual (also known as a *vector*) is equal to the number of rows of `bounds`.
That is, each row of `bounds` denotes the range ``[l, u]`` of a variable, where ``l`` and ``u`` denote 
the lower and upper bound, respectively. 

`TF` is the data type of fitness. The type of individuals (solution vectors) is the same as the type
of `bounds`.

The size of `senses` indicates the number of objectives.
If the supplied `senses` is a scalar, then a single-objective DE [`StdSoDE`](@ref) instance is built;
otherwise, a multi-objective [`StdMoDE`](@ref) is constructed. 
"""
function StdDE(::Type{TF}, nindividuals::Int, bounds::AM{<:AbstractFloat}, opt::DEOptions; 
        senses::ASoV{<:Integer}) where {TF<:Real}
    d = size(bounds, 1)  # number of dimensions
    # init a population randomly within bounds
    pop = zeros(eltype(bounds), d, nindividuals)
    @inbounds for i in 1:d
        l, u = bounds[i, :]
        @assert l <= u
        pop[i, :] .= l .+ rand(eltype(bounds), nindividuals) .* (u - l)
    end 
    # fitness set to default values zero
    if senses isa Real 
        fitness = zeros(TF, nindividuals)
    else
        fitness = zeros(TF, length(senses), nindividuals) 
    end
    return StdDE(pop, fitness, senses, bounds, opt)
end

"""
    StdDE(nindividuals::Int, bounds::AM{<:AbstractFloat}, opt::DEOptions; senses::ASoV{<:Integer})

Construct a standard DE whose vectors are of default type `Float64`.
See [`StdDE(::Type{TF}, nindividuals::Int, bounds::AM{<:AbstractFloat}, opt::DEOptions; 
        senses::ASoV{<:Integer}) where {TF<:Real}`](@ref) for details.
"""
StdDE(nindividuals::Int, bounds::AM{<:AbstractFloat}, opt::DEOptions; senses::ASoV{<:Integer}) =
    StdDE(Float64, nindividuals, bounds, opt; senses)

"""
    StdSoDE{TP<:AM{<:AF}, TF<:AV{<:Real}, TW<:Real, TO<:DEOptions}

Standard single-objective DE.
"""
const StdSoDE{TP<:AM{<:AF}, TF<:AV{<:Real}, TW<:Real, TO<:DEOptions} = StdDE{TP, TF, TW, TO}

"""
    StdMoDE{TP<:AM{<:AF}, TF<:AM{<:Real}, TW<:AV{<:Real}, TO<:DEOptions}

Standard multi-objective DE.
"""
const StdMoDE{TP<:AM{<:AF}, TF<:AM{<:Real}, TW<:AV{<:Real}, TO<:DEOptions} = StdDE{TP, TF, TW, TO}

"""
    best_individual(de::StdSoDE; all=false)

Get the best individual (or individuals) from the single-objective DE population.
See also [`best_individual(e::Evolver; all=false)`](@ref) for more details.
"""
function best_individual(de::StdSoDE; all=false)
    bf = best_fitness(de)
    if all
        ib = findall(isapprox(bf), fitness(de))
    else
        ib = findfirst(isapprox(bf), fitness(de))
    end
    return de.pop[:, ib]
end

"""
    best_fitness(de::StdSoDE)

Get the best fitness from the single-objective DE population.
"""
function best_fitness(de::StdSoDE)
    if de.senses > 0 # a scalar here
        return maximum(fitness(de))
    else
        return minimum(fitness(de))
    end
end

"""
    sort!(de::StdDES; rev::Bool=false)

Sort a single-objective DE population in place according to the fitness. The *best* individual appears
at first if `rev` is `false` (the default value).
"""
function Base.sort!(de::StdSoDE; rev::Bool=false)
    sort!(de._rand_indices; by = i -> de.fitness[i] * de.weights, rev=true)
    de.fitness .= de.fitness[de._rand_indices]
    de.pop .= de.pop[:, de._rand_indices]
end

# TODO Put the controller into evolve
Base.@kwdef struct StdSoController{T<:AF, C<:Union{Nothing, Function}}
    atol::T = 1e-20
    rtol::T = 1e-20
    successive::Int = 10
    ngen::Int = 500
    terminator::C = nothing
end


"""
    evolve!(de::StdDES, evaluator::Function, ngen; mut=mutate_rand_1, cx=crossover_binomial,
        bc=:bounce_back)

Standard single-objective differential evolution. 

Each individual's fitness is evaluated by the given `evaluator` and the evolution lasts for `ngen`
iterations. Keyword arguments:
- `mut=mutate_rand_1!`: mutation method
- `cx=crossover_binomial!`: crossover method, currently support
- `bc=:bounce_back`: boundary constraints handling method， currently support
    + :bounce_back
- `stats=[:best, :avg]`: fitness statistics to display in each generation. 
    Valid entries are `:best`, `:avg`(average), `:max`, `:min`, `:std`(standard deviation).
- `stats_fmt::FormatSpec=FormatSpec("10.5e")`: format specification of each statistic term for display.
    See [Formatting.jl](https://github.com/JuliaIO/Formatting.jl) for details.
- `history::Bool`: if `true`, record the convergence history including the max/min/avg/std fitness in 
    each interation and return the record as an array
"""
function evolve!(de::StdSoDE, evaluator::Function, ngen; mut=mutate_rand_1!, cx=crossover_binomial!,
    bc::Symbol=:bounce_back, stats=[:best, :avg], stats_fmt::FormatSpec=FormatSpec("<10.5e"),
    history::Bool=false)

    stats_func = Dict(:best => fs -> senses(de) > 0 ? maximum(de.fitness) : minimum(de.fitness),
        :avg=>mean, :std=>std, :max=>maximum, :min=>minimum)
    to_display = !isempty(stats)
    # display header 
    if to_display
        gen_width = 5
        stats_width = stats_fmt.width
        printfmt("{1:<$(gen_width)s}\t", "gen")
        for stat in stats
            printfmt("{1:<$(stats_width)s}\t", string(stat))
        end
        println()
    end
    # display each field 
    function report(g)
        printfmt("{1:<$(gen_width)d}\t", g)
        for stat in stats
            printfmt(stats_fmt, stats_func[stat](de.fitness))
            print("\t")
        end
        println()
    end

    n = nindividuals(de)
    s = senses(de)
    # initial evaluation
    de.fitness .= evaluator.(eachcol(de.pop))
    to_display && report(0)
    # main loop
    V = de.pop[:, 1]  # copy, to avoid frequent allocation
    U = de.pop[:, 1]
    # history logging
    if history
        records = zeros(ngen, 4)
    else
        records = zeros(0, 0)
    end
    for g in 1:ngen
        for i in 1:n
            mut(V, i, de) # donor
            # handle boundary constraints
            if bc == :bounce_back
                bounce_back!(V, view(de.pop, :, i), bounds(de))
            else
                error("Unimplemented")  # TODO
            end
            cx(U, i, de, V)  # trial
            obj = evaluator(U)
            if obj * s >= de.fitness[i] * s
                de.pop[:, i] .= U
                de.fitness[i] = obj
            end
        end
        to_display && report(g)
        if history
            # max/min/avg/std
            for (j, op) in enumerate((maximum, minimum, mean, std))
                records[g, j] = op(de.fitness)
            end
        end
        g += 1
    end
    return records
end
