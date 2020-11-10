#===================================================================================================
Abstract types and common interfaces.
===================================================================================================#
abstract type Evolver end
abstract type Options end
abstract type Controller end

# A normal (single-population) DE evolution type 
abstract type DEEvolver <: Evolver end

# Base type of options for a DE, particually, a single-population one. See also [`DEEvolver`](@ref).
abstract type DEOptions <: Options end
abstract type DEController <: Controller end

_throw_unimplemented_error() = error("Unimplemented")

# all subtypes of `DEEvolver` should have the following fields
pop(e::DEEvolver) = e.pop
fitness(e::DEEvolver) = e.fitness
senses(e::DEEvolver) = e.senses
bounds(e::DEEvolver) = e.bounds
options(e::DEEvolver) = e.options


"""
    nindividuals(e::DEEvolver) 

Number of individuals in the population.
"""
function nindividuals(e::DEEvolver) 
    return size(pop(e), 2)
end

"""
    ndims(de::DEEvolver)

Number of variables (dimension) of the DE.
"""
Base.ndims(de::DEEvolver) = size(pop(de), 1)

"""
    nobjectives(e::Evolver)

Number of objectives.
"""
function nobjectives(e::DEEvolver)
    f = fitness(e)
    if ndims(f) == 1
        return 1
    end
    return size(f, 1)
end

"""
    best_individual(e::Evolver; all=false)

Get individuals that have the best fitness. If `all` is `true`, then all such individuals are returned;
otherwise, only one of them is fetched.
"""
function best_individual(e::Evolver; all=false)
    _throw_unimplemented_error()
end

"""
    best_fitness(e::Evolver)

Get the best fitness among the individuals of the evolver.
"""
function best_fitness(e::Evolver)
    _throw_unimplemented_error()
end


