
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

function crossover_PMX(dst :: Vector{Int64}, src :: Vector{Int64}, idst :: Vector{Int64})
    # Build the inverse permutation for use as a lookup table.
    # Effectively working as a find-city but it does need to be kept up-to-date with each swap.
    invperm!(idst, dst)
    # Pick two points on the string
    ms_a, ms_b = rand(1:length(dst)), rand(1:length(dst))
    # These two points from a section on the string (which is assumed to be a ring)
    # eg. if b < a, it wraps around the end.
    if ms_b >= ms_a
        for i in ms_a:ms_b
            crossover_PMX_point(dst, src, idst, i)
        end
    else
        for i in ms_a:length(dst)
            crossover_PMX_point(dst, src, idst, i)
        end
        for i in 1:ms_b
            crossover_PMX_point(dst, src, idst, i)
        end
    end
    return dst
end

function crossover_CX(dst :: Vector{Int64}, src :: Vector{Int64}, idst :: Vector{Int64}, isrc :: Vector{Int64})
    # Create lookup tables for finding cycles.
    invperm!(idst, dst)
    invperm!(isrc, src)
    # Get a starting point, and store it, so we know when we have found a cycle!
    c = rand(1:length(dst))
    # Go through the cycle, copying over elements from src.
    i = c
    while isrc[dst[i]] != c
        # TODO: A variant that creates both offspring in one go?
        pi = i
        i = isrc[dst[i]]
        dst[pi] = src[pi]
    end
    dst[i] = src[i]
    
    return dst
end
