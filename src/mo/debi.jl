
using Distances: pairwise, Euclidean, euclidean
using Statistics: mean

# TODO: inplace
function cal_novelty(X::AM{<:AF})
    D = pairwise(Euclidean(), X, dims=2)
    return mean(D; dims=2)
end

function cal_novelty(x::AV{<:AF}, X::AM{<:AF})
    @assert length(x) == size(X, 1)
    return mean(y->euclidean(x, y), eachcol(X))
end

function debi!(de::StdMoDE, evaluator::Function, ngen; show_progress::Bool=true)::Vector{Vector{Int}}
    n = nindividuals(de)
    d = ndims(de)
    ss = senses(de)
    # initial evaluation
    de.fitness[1, :] .= evaluator.(eachcol(de.pop))
    de.fitness[2, :] .= vec(cal_novelty(de.pop))
    # de.fitness[2, :] .= 0
    # in DEMO, we will have an intermediate population of size [n, 2n]
    im_pop = zeros(eltype(de.pop), d, 2*n)
    im_fitness = zeros(eltype(de.fitness), nobjectives(de), 2*n)
    im_indices = collect(1:2*n)
    U = de.pop[:, 1] # copy to avoid allocation
    V = de.pop[:, 1]
    F = de.options.F
    Cr = de.options.Cr
    if show_progress
        pg = Progress(ngen, 0.1)
        ProgressMeter.ijulia_behavior(:clear)
    end
    fit = zeros(2)
    archive = de.pop
    for g = 1:ngen
        im_pop[:, 1:n] .= de.pop
        im_fitness[:, 1:n] .= de.fitness
        im_indices .= 1: 2*n
        cursor = n # the actual number of individuals in `im_pop`
        mean_obj = mean(@view im_fitness[1, :])
        for i = 1:n  # produce n offspring
            mutate_rand_1!(V, i, im_pop, F, @view im_indices[1:cursor]) # donor
            X = view(im_pop, :, i)  # target (base)
            bounce_back!(V, X, de.bounds)
            crossover_binomial!(U, X, V, Cr)  # trial
            fit[1] = evaluator(U)
            if fit[1] > mean_obj
                # continue
            end
            
            fit[2] = cal_novelty(U, view(im_pop, :, 1:cursor))
            # fit[2] = 0
            parent_fit = @view im_fitness[:, i]
            if dominate(fit, parent_fit, ss) # better, replace the parent
                im_pop[:, i] .= U
                im_fitness[:, i] .= fit
            elseif !dominate(parent_fit, fit, ss)  # equally good, add to the population
                cursor += 1
                im_pop[:, cursor] .= U
                im_fitness[:, cursor] .= fit
            end
        end
        # now we have a population of size `cursor` ∈ [n, 2n]
        # follow NSGA-II to reduce the population to size n
        # println("$g", ": $cursor")
        if cursor > n
            fronts = nondominated_sort(@view(im_fitness[:, 1:cursor]), ss)
            i_ind = 0
            i_front = 1
            # copy individuals in preferential fronts to de
            while i_ind + length(fronts[i_front]) <= n
                front = fronts[i_front]
                l = length(front)
                de.pop[:, i_ind+1:i_ind+l] .= @view im_pop[:, front]
                de.fitness[:, i_ind+1:i_ind+l] .= @view im_fitness[:, front]
                i_ind += l
                i_front += 1
            end
            # if there are still vacancies
            # we need only a portion of the next front accroding to crowding distance
            if i_ind < n
                front = fronts[i_front]
                l = length(front)
                im_indices[1:l] .= 1:l
                cds_i_front = assign_crowding_distance(view(im_fitness, :, front), [@view im_indices[1:l]])
                # arg-sort of the front in decending order w.r.t crowding distance
                ix = @view im_indices[1:l] # just a buffer
                sortperm!(ix, cds_i_front; rev=true)
                # pick the ones with larger crowding distance
                for i = 1:n-i_ind
                    picked = front[ix[i]]
                    de.pop[:, i_ind+i] .= @view im_pop[:, picked]
                    de.fitness[:, i_ind+i] .= @view im_fitness[:, picked]
                end
            end
        else # we have coincidentally n individuals
            front = Vector{Int}[]
            de.pop .= @view im_pop[:, 1:n]
            de.fitness .= @view im_fitness[:, 1:n]
        end

        if show_progress
            ProgressMeter.next!(pg; showvalues=[(:gen, "$g/$ngen")])
        end
    end
    # finally, perform nondominated sorting again to obtain fronts
    return nondominated_sort(de.fitness, ss)
end
