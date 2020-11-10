

"""
    crossover_binomial(i::Integer, de::StdDE{<:AM{T}}, V::AV{T}) where T<:AF
    crossover_binomial!(U::AV{T}, i::Integer, de::StdDE{<:AM{T}}, V::AV{T}) where T<:AF

Standard binomial crossover that recombines the donor vector `V` with the `i`-th target vector in 
the population of `de` to yield a donor vector.
"""
function crossover_binomial(i::Integer, de::StdDE{<:AM{T}}, V::AV{T}) where T<:AF
    U = de.pop[:, 1] # copy
    crossover_binomial!(U, i, de, V)
    return U
end

function crossover_binomial!(U::AV{T}, i::Integer, de::StdDE{<:AM{T}}, V::AV{T}) where T<:AF
    crossover_binomial!(U, view(de.pop, :, i), V, de.options.Cr)
end

function crossover_binomial!(U::AV{T}, X::AV{T}, V::AV{T}, Cr::AF) where T<:AF
    @assert length(U) == length(X) == length(V)
    d = length(U)
    j_rand = rand(1:d)
    @inbounds for j = 1:d
        if rand(typeof(Cr)) <= Cr || j == j_rand
            U[j] = V[j]
        else
            U[j] = X[j]
        end
    end
end
