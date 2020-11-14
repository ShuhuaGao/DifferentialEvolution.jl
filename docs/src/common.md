# Common interfaces of DE

## Core DE
The core data structure of differential evolution is a `struct` that inherits the abstract type
[`Evolver`](@ref) or, more typically, the [`DEEvolver`](@ref), which itself inherits `Evolver` and 
represents a single-population DE. 
```@docs
Evolver
Options
DEEvolver
DEOptions
```

Each concrete type of `DEEvolver` should support the five default interfaces listed below. If the five corresponding fields are defined, then nothing extra needs to be done.
```julia
pop(e::DEEvolver) = e.pop
fitness(e::DEEvolver) = e.fitness
senses(e::DEEvolver) = e.senses
bounds(e::DEEvolver) = e.bounds
options(e::DEEvolver) = e.options
```
Thus, the five fundamental fields listed above can be retrieved directly or using a method with the same name of the field.

```@docs
pop
fitness
senses
bounds
options
```


## DE options
Options, i.e., control parameters, of DE are encapsulated in a `struct` as a subtype of [`DEOptions`](@ref) (for the most common single-population DE) or a subtype of [`Options`](@ref). The options of a given DE is accessed by the method [`options`](@ref).

## Other utility interfaces
The following three methods have default implementation for [`DEEvolver`](@ref).
```@docs
nindividuals(e::DEEvolver) 
ndims(de::DEEvolver)
nobjectives(e::DEEvolver)
```

No default implementation is provided for the two methods below, since a correct implementation depends on details of the concrete DE type.
```@docs
best_individual(e::Evolver; all=false)
best_fitness(e::Evolver)
```
It is highly recommended that a concrete DE type implements the two *best* methods as outlined above.

Currently, `DifferentialEvolution.jl` provides a *standard* implementation with the types [`StdDE`](@ref) and [`StdOptions`](@ref).