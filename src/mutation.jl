#===================================================================================================
Mutation operators in differential evolution.
===================================================================================================#

"""
    mutate_rand_1(i::Integer, de::StdDE)

DE/rand/1 mutation, which produces a donor vector to crossover with the `i`-th target vector (parent).
"""
function mutate_rand_1(i::Integer, de::StdDE)
    V = de.pop[:, 1]  # copy 
    mutate_rand_1!(V, i, de)
    return V
end

"""
    mutate_rand_1!(V::AV{T}, i::Integer, de::StdDE{<:AM{T}}) where T<:AF

DE/rand/1 mutation in place, which produces a donor vector to crossover with the `i`-th target vector (parent).
"""
function mutate_rand_1!(V::AV{T}, i::Integer, de::StdDE{<:AM{T}}) where T<:AF
    mutate_rand_1!(V, i, de.pop, de.options.F, de._rand_indices)
end

 #TODO current

function mutate_rand_1!(V::AV{T}, i::Integer, pop::AM{T}, F::AF, indices::AV{<:Integer}) where T<:AF
    shuffle!(indices)
    # choose three vectors different from `i` for differential operation    
    chosen = MVector{3, eltype(indices)}(0, 0, 0)
    count = 0
    for idx in indices
        if idx != i
            count += 1
            chosen[count] = idx
            if count == 3
                break
            end
        end
    end
    r1, r2, r3 = chosen
    @views V .= pop[:, r1] .+ F.*(pop[:, r2] .- pop[:, r3])
    
end