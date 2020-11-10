#===================================================================================================
Infrastructure of multiobjective evolution. 
===================================================================================================#

"""
    dominate(f1::AV{T}, f2::AV{T}, [senses::AV{<:Integer}]) where {T<:Real}

Test whether `f1` dominates `f2` w.r.t the given `senses`. If the argument `senses` is not provided, 
then maximization is assumed.
"""
function dominate(f1::AV{T}, f2::AV{T}, senses::AV{<:Integer}) where {T<:Real}
    @assert length(f1) == length(f2) && length(f2) == length(senses)
    one_better = false
    m = length(f1)
    @inbounds for i = 1:m
        fw1 = f1[i] * senses[i] 
        fw2 = f2[i] * senses[i]
        if fw1 < fw2 
            return false
        elseif fw1 > fw2
            one_better = true
        end
    end
    return one_better
end

function dominate(f1::AbstractVector{T}, f2::AbstractVector{T}) where {T<:Real}
    @assert length(f1) == length(f2)
    one_better = false
    m = length(f1)
    @inbounds for i = 1:m
        fw1 = f1[i]
        fw2 = f2[i]
        if fw1 < fw2 
            return false
        elseif fw1 > fw2
            one_better = true
        end
    end
    return one_better
end