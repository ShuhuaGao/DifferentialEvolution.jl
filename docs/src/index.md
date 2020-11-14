# Documentation of DifferentialEvolution.jl

## Type aliases
To avoid writing lengthy type names, `DifferentialEvolution.jl` defines the following abbreviations.
```julia
const AbstractScalarOrVec{T} = Union{T, AbstractVector{T}} where T<:Real
const AV = AbstractVector
const AM = AbstractMatrix
const AVoM = AbstractVecOrMat
const ASoV = AbstractScalarOrVec
const AF = AbstractFloat
```

## Contents

[Common interfaces of DE](@ref)

[Mutation](@ref)

[Crossover](@ref)

[Standard DE](@ref)