#===================================================================================================
Constraint handling differential evolution
===================================================================================================#
"""
    bounce_back!(t::AV{T}, b::AV{T}, bounds::AM{<:Real}) where T<:AF

Relocates randomly the parameter of the trial vector `t` in between the bound it exceeded and the 
corresponding parameter from the base vector `b`.
"""
function bounce_back!(t::AV{T}, b::AV{T}, bounds::AM{<:Real}) where T<:AF
    @assert length(t) == length(b) == size(bounds, 1)
    @inbounds for i in eachindex(t)
        lower, upper = bounds[i, 1], bounds[i, 2]
        if t[i] < lower  # smaller than the lower bound
            t[i] = lower + rand(T) * (b[i] - lower)
        elseif t[i] > upper # exceed the upper bound
            t[i] = upper - rand(T) * (upper - b[i])
        end
    end
end
