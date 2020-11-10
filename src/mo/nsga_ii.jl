#===================================================================================================
NSGA-II

[1] Deb, Kalyanmoy, et al. "A fast and elitist multiobjective genetic algorithm: NSGA-II." 
IEEE transactions on evolutionary computation 6.2 (2002): 182-197.
===================================================================================================#

"""
    nondominated_sort(fitness::AbstractMatrix{<:Real}, senses::AbstractVector{<:Integer})::Vector{Vector{Int}}

Fast nondominated sorting on page 3 of [1]. The fronts are returned, where the first one is the Pareto
optimal front, and the i-th front contains the indices of individuals in that front.  
"""
function nondominated_sort(fitness::AbstractMatrix{<:Real}, senses::AbstractVector{<:Integer})::Vector{Vector{Int}}
    n = size(fitness, 2) # number of individuals
    d_counter = zeros(Int, n)
    dominated = [Int[] for _ = 1:n]
    fronts = Vector{Int}[]
    @inbounds for ip = 1:n, iq = ip+1:n
        p = view(fitness, :, ip)
        q = view(fitness, :, iq)
        if dominate(p, q, senses)
            push!(dominated[ip], iq)
            d_counter[iq] += 1
        elseif dominate(q, p, senses)
            push!(dominated[iq], ip)
            d_counter[ip] += 1
        end
    end
    # the first front
    push!(fronts, Int[])
    f1 = fronts[1]
    @inbounds for i = 1:n
        if d_counter[i] == 0
            push!(f1, i)
        end
    end
    n_in_fronts = length(f1)
    # the other fronts
    current = 1
    while true
        new_front = Int[]
        sizehint!(new_front, n - n_in_fronts)  
        for ip in fronts[current]
            for iq in dominated[ip]
                d_counter[iq] -= 1
                if d_counter[iq] == 0
                    push!(new_front, iq)
                end
            end
        end
        if isempty(new_front)
            break
        else
            push!(fronts, new_front)
            n_in_fronts += length(new_front)
        end
        current += 1
    end
    return fronts
end


"""
    get_nondomination_rank(fronts::AbstractVector{<:AbstractVector{<:Integer}})::Vector{Int}

Given the sorted front, determine which front each individual resides (starting from 1).
The value at the i-th position of the returned array is the rank of the i-th individual in the 
population.
"""
function get_nondomination_rank(fronts::AbstractVector{<:AbstractVector{<:Integer}})::Vector{Int}
    n = sum(length, fronts)
    ranks = zeros(Int, n)
    for (r, front) in enumerate(fronts)
        ranks[front] .= r
    end
    return ranks
end


"""
    assign_crowding_distance(fitness::AbstractMatrix{<:Real}, 
        fronts::AbstractVector{<:AbstractVector{<:Integer}})::Vector{Real}

Compute the crowding distance for each individual in the whole population characterized by its 
`fitness` and Pareto front in  `fronts`. Refer to the algorithm on the bottom of page 4 of Ref. [1].

A vector of crowding distances `cds` is returned, in which `cds[i]` denotes the crowding distance of 
the individual `i`.
"""
function assign_crowding_distance(fitness::AbstractMatrix{<:Real}, 
    fronts::AbstractVector{<:AbstractVector{<:Integer}})::Vector{Float64}
    m, n = size(fitness)
    cds = zeros(n)  # crowding distance
    fronts = deepcopy(fronts)
    # handle each front separately
    for front in fronts # front contains the index of an individual in the population
        if length(front) <= 2   
            cds[front] .= Inf
            continue
        end
        for io = 1:m  # each objective
            sort!(front; by = idx -> fitness[io, idx])
            scale = fitness[io, front[end]] - fitness[io, front[1]]
            # all individuals in this front have the identical `io`-th objective value, ignore it
            if scale > 0  
                cds[front[1]] = Inf
                cds[front[end]] = Inf
                for i = 2:length(front)-1
                    cds[front[i]] += (fitness[io, front[i+1]] - fitness[io, front[i-1]]) / scale
                end
            end
        end
    end
    return cds    
end


"""
   crowded_compare(i_rank, i_dist, j_rank, j_dist)::Bool
   
Crowed-comparison operator. Return -1 if we prefer the individual `i`; otherwise, 1 means
the individual `j` is preferable. If both individuals `i` and `j` are equally good, 0 is returned.
"""
function crowded_compare(i_rank, i_dist, j_rank, j_dist)::Int
    if i_rank < j_rank
        return -1
    elseif i_rank > j_rank
        return 1
    end
    # now i and j belong to the same Pareto front: prefer the less crowded one
    if i_dist > j_dist
        return -1
    elseif i_dist < j_dist
        return 1
    end
    return 0
end


