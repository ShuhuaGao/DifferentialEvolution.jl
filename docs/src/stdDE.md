# Standard DE
The standard single-objective DE types and evolution algortithms. 

## Types

```@docs
StdOptions
StdDE
StdSoDE
StdMoDE
```

## Constructor
The three standard DE types [`StdDE`](@ref), [`StdSoDE`](@ref), and [`StdMoDE`](@ref) can all be built with their inner constructors, but a more convenient way is to employ the dedicated outer constructor. Please check Julia documentation on [Constructors](https://docs.julialang.org/en/v1/manual/constructors/) for more details.

```@meta
CurrentModule = DifferentialEvolution
```

```@docs
StdDE(::Type{TF}, nindividuals::Int, bounds::AM{<:AbstractFloat}, opt::DEOptions; 
        senses::ASoV{<:Integer}) where {TF<:Real}

StdDE(nindividuals::Int, bounds::AM{<:AbstractFloat}, opt::DEOptions; senses::ASoV{<:Integer})
```
## Evolution
The basic single-objective DE evolution algorithm is implemented.
```@docs
evolve!
```

## Ranking
```@docs
best_individual(de::StdSoDE; all=false)
best_fitness(de::StdSoDE)
Base.sort!(de::StdSoDE; rev::Bool=false)
```