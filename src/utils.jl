# abbreviations of types with long names
const AbstractScalarOrVec{T} = Union{T, AbstractVector{T}} where T<:Real
const AV = AbstractVector
const AM = AbstractMatrix
const AVoM = AbstractVecOrMat
const ASoV = AbstractScalarOrVec
const AF = AbstractFloat