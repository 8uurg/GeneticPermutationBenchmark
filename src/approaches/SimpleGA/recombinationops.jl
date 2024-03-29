## Crossover Implementations.
function invperm!(dst :: Vector{Int64}, src :: Vector{Int64})
    @inbounds for (i, j) in enumerate(src)
        dst[j] = i
    end
    return #dst
end

function crossover_PMX_point(dst :: Vector{Int64}, src :: Vector{Int64}, idst :: Vector{Int64}, i :: Int64)
    # What numbers are we swapping
    a = src[i]
    b = dst[i]
    # Perform the swap
    dst[i], dst[idst[a]] = dst[idst[a]], b
    # Swap the numbers in the lookup table as well.
    idst[a], idst[b] = idst[b], idst[a]
    return
end

function crossover_PMX!(dst1 :: Vector{Int64}, dst2 :: Vector{Int64}, parent1 :: Vector{Int64}, parent2 :: Vector{Int64},
        idst1 :: Vector{Int64}, idst2 :: Vector{Int64}, rng :: MersenneTwister, ring :: Bool = true)
    n = length(parent1)
    #
    copy!(dst1, parent1)
    copy!(dst2, parent2)
    # Build the inverse permutation for use as a lookup table.
    # Effectively working as a find-city but it does need to be kept up-to-date with each swap.
    invperm!(idst1, dst1)
    invperm!(idst2, dst2)
    # Pick two points on the string
    ms_a, ms_b = rand(rng, 1:n), rand(rng, 1:n)
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
    return #dst1, dst2
end

function crossover_CX!(dst1 :: Vector{Int64}, dst2 :: Vector{Int64}, parent1 :: Vector{Int64}, parent2 :: Vector{Int64},
        idst2 :: Vector{Int64}, rng :: MersenneTwister)
    n = length(parent1)
    # Offspring start of as each of the parents.
    copy!(dst1, parent1)
    copy!(dst2, parent2)
    # Create lookup table for finding cycles.
    invperm!(idst2, dst2)
    # Get a starting point, and store it, so we know when we have found a cycle!
    c = rand(rng, 1:n)
    # Go through the cycle, copying over elements from src.
    i = c
    while idst2[dst1[i]] != c
        pi = i
        i = idst2[dst1[i]]
        dst1[pi], dst2[pi] = dst2[pi], dst1[pi]
    end
    dst1[i], dst2[i] = dst2[i], dst1[i]

    return #dst1, dst2
end

function crossover_OX!(dst1 :: Vector{Int64}, dst2 :: Vector{Int64}, parent1 :: Vector{Int64}, parent2 :: Vector{Int64},
        seen_dst1 :: BitVector, seen_dst2 :: BitVector, rng :: MersenneTwister, ring :: Bool = true)
    n = length(parent1)
    # Select crossover points
    ms_a, ms_b = rand(rng, 1:n), rand(rng, 1:n)
    # If not acting like a ring, reorder the points.
    if !ring && ms_a > ms_b
        ms_a, ms_b = ms_b, ms_a
    end
    # Keep track of which items we have seen
    fill!(seen_dst1, false)
    fill!(seen_dst2, false)
    # Case distinction between wrapping around and
    # non-wraparound.
    if ms_b >= ms_a
        # Copy over elements in window
        for i in ms_a:ms_b
            dst1[i] = parent1[i]
            seen_dst1[parent1[i]] = true
            dst2[i] = parent2[i]
            seen_dst2[parent2[i]] = true
        end
        # From other parent, copy over elements
        # for positions outside of window.
        p1 = 1
        p2 = 1
        for i in 1:ms_a-1
            # Ignore seen elements.
            while seen_dst2[parent1[p1]]
                p1 += 1
            end
            while seen_dst1[parent2[p2]]
                p2 += 1
            end
            dst2[i] = parent1[p1]
            dst1[i] = parent2[p2]
            p1 += 1
            p2 += 1
        end
        for i in ms_b+1:n
            # Ignore seen elements.
            while seen_dst2[parent1[p1]]
                p1 += 1
            end
            while seen_dst1[parent2[p2]]
                p2 += 1
            end
            dst2[i] = parent1[p1]
            dst1[i] = parent2[p2]
            p1 += 1
            p2 += 1
        end
    else
        # Copy over elements in window
        for i in ms_a:n
            dst1[i] = parent1[i]
            seen_dst1[parent1[i]] = true
            dst2[i] = parent2[i]
            seen_dst2[parent2[i]] = true
        end
        for i in 1:ms_b
            dst1[i] = parent1[i]
            seen_dst1[parent1[i]] = true
            dst2[i] = parent2[i]
            seen_dst2[parent2[i]] = true
        end
        p1 = 1
        p2 = 1
        for i in ms_b+1:ms_a-1
            # Ignore seen elements.
            while seen_dst2[parent1[p1]]
                p1 += 1
            end
            while seen_dst1[parent2[p2]]
                p2 += 1
            end
            dst2[i] = parent1[p1]
            dst1[i] = parent2[p2]
            p1 += 1
            p2 += 1
        end
    end
    return #dst1, dst2
end

function crossover_ER!(dst :: Vector{Int64}, parent1 :: Vector{Int64}, parent2 :: Vector{Int64}, adj :: Matrix{Int64}, cnt :: Vector{Int64}, remaining :: Vector{Int64}, remaining_lut :: Vector{Int64}, rng :: MersenneTwister, ring :: Bool = true)
    #
    n = length(parent1)
    # Prepare remaining list and lut.
    @inbounds for i in 1:n
        remaining[i] = i
        remaining_lut[i] = i
    end
    remaining_c = n
    # Create adjacency matrix, 2 parents = maximum of 4 adjacent elements.
    # Predecessor and successor in a single solution are always unique. (No need to check!)
    fill!(cnt, 0)
    fill!(adj, 0)
    @inbounds for i in 1:n
        x = parent1[i]
        for parentx in 1:2
            if parentx == 1
                p = parent1
            else
                p = parent2
            end
            x = p[i]
            pred = p[mod((i - 1 - 1), n) + 1]
            succ = p[mod((i + 1 - 1), n) + 1]
            if cnt[x] != 0
                # Duplicate checking with found adjacent elements.
                if pred != adj[x, 1] && pred != adj[x, 2] && (ring || i != 1)
                    cnt[x] += 1
                    adj[x, cnt[x]] = pred
                end
                if succ != adj[x, 1] && succ != adj[x, 2] && (ring || i != n)
                    cnt[x] += 1
                    adj[x, cnt[x]] = succ
                end
            else
                # Easier case. No chance of duplicates.
                # though we will need to account for ring behaviour.
                if ring || i != 1
                    cnt[x] += 1
                    adj[x, cnt[x]] = pred
                end
                if ring || i != n
                    cnt[x] += 1
                    adj[x, cnt[x]] = succ
                end
            end
        end
    end

    # Construct a solution. Starting from the starting point of one of the two parents
    # selected at random.
    @inbounds c = rand(rng, (parent1[1], parent2[1]))
    @inbounds dst[1] = c
    # Update adjacency matrix of each neighbour
    @inbounds for ni in 1:cnt[c]
        # Neighbour `n`
        nb = adj[c, ni]
        # Find index of ci.
        ci = 0
        for b in 1:cnt[nb]
            if adj[nb, b] == c
                ci = b
                break
            end
        end
        # Update adjacency and count
        if ci != 0 && cnt[nb] != 0
            adj[nb, ci] = adj[nb, cnt[nb]]
            cnt[nb] -= 1
        end
    end
    # Mark as selected
    @inbounds cnt[c] = -1
    # Move the last unselected item
    # to the position of the selected element
    # and 'mark' the last spot as selected.
    @inbounds remaining[remaining_lut[c]] = remaining[remaining_c]
    @inbounds remaining_lut[remaining[remaining_c]] = remaining_lut[c]
    remaining_c -= 1
    # Do the same thing for the remaining elements
    @inbounds for i in 2:n
        # Find the element with smallest count
        # and randomly pick from the set of the elements with smallest count.
        nxt = 0
        nxtc = 5
        w = 1
        for a in 1:cnt[c]
            p = adj[c, a]
            if cnt[p] < nxtc
                nxt = p
                nxtc = cnt[p]
                w = 1
            elseif cnt[p] == nxtc
                w += 1
                # Current item is the result of a sample over w-1 items.
                # Probability the next one should be selected is 1.0/w
                if rand(rng, Float64) < 1.0/convert(Float64, w)
                    nxt = p
                end
            end
        end

        # No neighbouring elements...
        # Select an unselected one at random instead.
        if nxt == 0
            w = rand(rng, 1:remaining_c)
            nxt = remaining[w]
        end
        # Next element has been chosen!
        c = nxt
        dst[i] = c
        remaining[remaining_lut[c]] = remaining[remaining_c]
        remaining_lut[remaining[remaining_c]] = remaining_lut[c]
        remaining_c -= 1

        # Update adjacency matrix of each neighbour
        for ni in 1:cnt[c]
            # Neighbour `n`
            nb = adj[c, ni]
            # Find index of ci.
            ci = 0
            for b in 1:cnt[nb]
                if adj[nb, b] == c
                    ci = b
                    break
                end
            end
            # Update adjacency and count
            if ci != 0 && cnt[nb] != 0
                adj[nb, ci] = adj[nb, cnt[nb]]
                cnt[nb] -= 1
            end
        end
        # Mark as selected.
        cnt[c] = -1
    end
    return #dst
end

abstract type CrossoverOperator end

struct PMX <: CrossoverOperator
    idst1 :: Vector{Int64}
    idst2 :: Vector{Int64}

    function PMX(n :: Int64)
        new(collect(1:n), collect(1:n))
    end
end

function crossover!(operator :: PMX,
                    offspring1 :: Vector{Int64},
                    offspring2 :: Vector{Int64},
                    parent1 :: Vector{Int64},
                    parent2 :: Vector{Int64},
                    rng :: MersenneTwister)
    #
    crossover_PMX!(offspring1, offspring2, parent1, parent2, operator.idst1, operator.idst2, rng, true)
end

struct CX <: CrossoverOperator
    idst2 :: Vector{Int64}

    function CX(n :: Int64)
        new(collect(1:n))
    end
end

function crossover!(operator :: CX,
                    offspring1 :: Vector{Int64},
                    offspring2 :: Vector{Int64},
                    parent1 :: Vector{Int64},
                    parent2 :: Vector{Int64},
                    rng :: MersenneTwister)
    #
    crossover_CX!(offspring1, offspring2, parent1, parent2, operator.idst2, rng)
end

struct OX <: CrossoverOperator
    seen_dst1 :: BitVector
    seen_dst2 :: BitVector

    function OX(n :: Int64)
        new(falses(n), falses(n))
    end
end

function crossover!(operator :: OX,
                    offspring1 :: Vector{Int64},
                    offspring2 :: Vector{Int64},
                    parent1 :: Vector{Int64},
                    parent2 :: Vector{Int64},
                    rng :: MersenneTwister)
    #
    crossover_OX!(offspring1, offspring2, parent1, parent2, operator.seen_dst1, operator.seen_dst2, rng, true)
end

struct ER <: CrossoverOperator
    adj :: Matrix{Int64}
    cnt :: Vector{Int64}
    remaining :: Vector{Int64}
    remaining_lut :: Vector{Int64}

    function ER(n :: Int64)
        new(zeros(Int64, (n, 4)), zeros(Int64, n), collect(1:n), collect(1:n))
    end
end

function crossover!(operator :: ER,
                    offspring1 :: Vector{Int64},
                    offspring2 :: Vector{Int64},
                    parent1 :: Vector{Int64},
                    parent2 :: Vector{Int64},
                    rng :: MersenneTwister)
    # Two crossovers, as Edge Recombination Crossover does not really have
    # have an obvious way to obtain two 'opposite' offspring
    crossover_ER!(offspring1, parent1, parent2, operator.adj, operator.cnt, operator.remaining, operator.remaining_lut, rng, true)
    crossover_ER!(offspring2, parent1, parent2, operator.adj, operator.cnt, operator.remaining, operator.remaining_lut, rng, true)
    return #(offspring1, offspring2)
end
