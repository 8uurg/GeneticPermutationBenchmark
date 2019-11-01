# This file contains utilities for working with Knjazew's Deceptive Permutation Function
# such as Storage (Struct), Evaluation and creating a Black-Box function. 
# !!!note The Knjazew's Deceptive Permutation Function is a **maximization** problem.

"A struct representing an instance of Knjazew's Deceptive Permutation Function"
struct KDPInstance
    k :: Int64
    n :: Int64 
end

"""
    distance(block :: Vector{Int64}) :: Int64
    distance(block :: Vector{Int64}, lo :: Vector{Int64} = ones(Int64, length(block))) :: Int64

Compute the minimum amount of integers to move in order to obtain the 'sorted' sequence.
Or alternatively: the longest ascending subsequence.
"""
function perm_distance(block, start :: Int64 = 1, len :: Int64 = length(block), lo :: Vector{Int64} = ones(Int64, len)) :: Int64
    # Each integer is atleast part of a ascending subsequence of size 1
    # with itself!
    fill!(lo, one(Int64))
    #
    dist :: Int64 = 1

    # For each i, find the integer that should appear before it
    # with the current largest ascending subsequence.
    @inbounds for i in start:start+len-1, j in i-1:-1:start
        if block[j-start+1] < block[i-start+1] && lo[j-start+1] >= lo[i-start+1]
            # The resulting sequence is one larger
            lo[i-start+1] = lo[j-start+1] + 1
            # Update maximum
            if lo[i-start+1] > dist
                dist = lo[i-start+1]
            end
        end
    end
    return length(block) - dist
end

"""
    evaluate_kdp(instance :: KDPInstance, permutation :: Vector{Int64})

Evaluate `assignment` against KDPInstance `instance`, and 
return its objective value.
"""
function evaluate_kdp(instance :: KDPInstance, 
        permutation :: Vector{Int64},
        lo :: Vector{Int64} = ones(Int64, instance.k),
        ip :: Vector{Int64} = zeros(Int64, instance.n*instance.k)) :: Float64
    # Sanity check: permutation should be of same length.
    @assert length(permutation) == instance.n * instance.k
    # Initially score is zero.
    score :: Float64 = zero(Float64)
    # Gather blocks together by taking the inverse permutation.
    # [1 .. k] ..(n times).. [1 .. k]
    # Note: if this is not layout you want, apply an additional transformation during this step
    for (i, p) in enumerate(permutation)
        ip[p] = i 
    end
    # Additively compose the fitness the blocks.
    # Note that perm_distance(a) == perm_distance(invperm(a)), hence perm_distance can be applied to the slice directly.
    for i in 1:instance.n
        dist = perm_distance(ip, (i-1)*instance.k+1, instance.k, lo)
        # The function is deceptive, this needs to be maximized.
        deceptive_dist = ifelse(dist > 0, convert(Float64, dist) / convert(Float64, instance.k), 1.0)
        score += deceptive_dist
    end

    return score
end

"""
    bb_wrap_kdp(instance :: KDPInstance)

Wrap an instance of the Knjazew's Deceptive Permutation Problem and return a black box
that takes an assignment and returns its fitness value on the given instance.
"""
function bb_wrap_kdp(instance :: KDPInstance)
    # Preallocate memory
    lo :: Vector{Int64} = ones(Int64, instance.k)
    ip :: Vector{Int64} = zeros(Int64, instance.n*instance.k)
    
    function evaluate(assignment :: Vector{Int64}) :: Float64
        return evaluate_kdp(instance, assignment, lo, ip)
    end
    return evaluate
end
