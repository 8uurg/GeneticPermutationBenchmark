
function invperm!(dst :: Vector{Int64}, src :: Vector{Int64})
    @inbounds for (i, j) in enumerate(src)
        dst[j] = i
    end
    dst
end

function crossover_PMX_point(dst :: Vector{Int64}, src :: Vector{Int64}, idst :: Vector{Int64}, i :: Int64)
    # What numbers are we swapping
    a = src[i]
    b = dst[i] 
    # Perform the swap
    dst[i], dst[idst[a]] = dst[idst[a]], b
    # Swap the numbers in the lookup table as well.
    idst[a], idst[b] = idst[b], idst[a]
end

function crossover_PMX_both!(dst1 :: Vector{Int64}, dst2 :: Vector{Int64}, parent1 :: Vector{Int64}, parent2 :: Vector{Int64}, 
        idst1 :: Vector{Int64}, idst2 :: Vector{Int64}; ring :: Bool = true)
    n = length(parent1)
    #
    copy!(dst1, parent1)
    copy!(dst2, parent2)
    # Build the inverse permutation for use as a lookup table.
    # Effectively working as a find-city but it does need to be kept up-to-date with each swap.
    invperm!(idst1, dst1)
    invperm!(idst2, dst2)
    # Pick two points on the string
    ms_a, ms_b = rand(1:n), rand(1:n)
    # If not acting like a ring, do a swap if neccesary.
    if !ring && ms_b < ms_a
        ms_a, ms_b = ms_b, ms_a
    end
    # These two points from a section on the string (which is assumed to be a ring)
    # eg. if b < a, it wraps around the end.
    if ms_b >= ms_a
        for i in ms_a:ms_b
            crossover_PMX_point(dst1, parent2, idst1, i)
            crossover_PMX_point(dst2, parent1, idst2, i)
        end
    else
        for i in ms_a:n
            crossover_PMX_point(dst1, parent2, idst1, i)
            crossover_PMX_point(dst2, parent1, idst2, i)
        end
        for i in 1:ms_b
            crossover_PMX_point(dst1, parent2, idst1, i)
            crossover_PMX_point(dst2, parent1, idst2, i)
        end
    end
    return dst1, dst2
end

function crossover_CX!(dst1 :: Vector{Int64}, dst2 :: Vector{Int64}, parent1 :: Vector{Int64}, parent2 :: Vector{Int64},
        idst2 :: Vector{Int64})
    n = length(parent1)
    # Offspring start of as each of the parents.
    copy!(dst1, parent1)
    copy!(dst2, parent2)
    # Create lookup table for finding cycles.
    invperm!(idst2, dst2)
    # Get a starting point, and store it, so we know when we have found a cycle!
    c = rand(1:n)
    # Go through the cycle, copying over elements from src.
    i = c
    while idst2[dst1[i]] != c
        pi = i
        i = idst2[dst1[i]]
        dst1[pi], dst2[pi] = dst2[pi], dst1[pi]
    end
    dst1[i], dst2[i] = dst2[i], dst1[i]
    
    return dst1, dst2
end
