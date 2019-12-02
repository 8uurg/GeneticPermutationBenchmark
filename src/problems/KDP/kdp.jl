# This file contains utilities for working with Knjazew's Deceptive Permutation Function
# such as Storage (Struct), Evaluation and creating a Black-Box function. 
# !!!note The Knjazew's Deceptive Permutation Function is a **maximization** problem.

# TODO: Make this work with instances for which n != b * k
# Eg. instances with partial blocks.

"A struct representing an instance of Knjazew's Deceptive Permutation Function"
struct KDPInstance
    k :: Int64 # size of a block
    b :: Int64 # number of blocks
    n :: Int64 # Number of elements
    mapping :: Vector{Int64}

    function KDPInstance(k :: Int64, b :: Int64)
        new(k, b, k*b, gen_block_mapping(k, b))
    end

    function KDPInstance(k :: Int64, n :: Int64, mapping :: Vector{Int64})
        new(k, n, k*b, mapping)
    end
end

"""
    gen_block_mapping(k :: Int64, b :: Int64)

Generate the mapping that maps neighbouring elements in
blocks to neighbouring elements in the actual vector.
Such that the same encoding is used as in
 "Design and Application of Iterated Density Evolutionary Algorithms"
  by P. A. N. Bosman

Here `k` is the number of elements in a block and `b` the number of blocks. 
"""
function gen_block_mapping(k :: Int64, b :: Int64)
    mapping = collect(1:k*b)
    i = 1
    for block in 1:b
        for element in 1:k
            mapping[block + b*(element-1)] = i 
            i += 1
        end
    end
    return mapping
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
    for i in start:start+len-1, j in i-1:-1:start
        # Stop early upon encountering a -1.
        if block[i] == -1
            return i - start - dist
        end
        if block[j] < block[i] && lo[j-start+1] >= lo[i-start+1]
            # The resulting sequence is one larger
            lo[i-start+1] = lo[j-start+1] + 1
            # Update maximum
            if lo[i-start+1] > dist
                dist = lo[i-start+1]
            end
        end
    end
    return len - dist
end

"""
    evaluate_kdp(instance :: KDPInstance, permutation :: Vector{Int64})

Evaluate `assignment` against KDPInstance `instance`, and 
return its objective value.
"""
function evaluate_kdp(instance :: KDPInstance, 
        permutation :: Vector{Int64},
        lo :: Vector{Int64} = ones(Int64, instance.k),
        ip :: Vector{Int64} = zeros(Int64, instance.n)) :: Tuple{Float64, Int64}
    # Sanity check: permutation should be of same length.
    @assert length(permutation) == instance.n
    # Initially score is zero.
    score :: Float64 = zero(Float64)
    blocks :: Int64 = 0
    # Gather blocks together by taking the inverse permutation.
    # [1 .. k] ..(n times).. [1 .. k]
    # Note: if this is not layout you want, apply an additional transformation during this step
    fill!(ip, -1)
    for (i, p) in enumerate(permutation)
        ip[instance.mapping[p]] = i
    end
    # Additively compose the fitness the blocks.
    # Note that perm_distance(a) == perm_distance(invperm(a)), hence perm_distance can be applied to the slice directly.
    for i in 1:instance.b
        dist = perm_distance(ip, (i-1)*instance.k+1, instance.k, lo)
        # The function is deceptive, this needs to be maximized.
        deceptive_dist = ifelse(dist > 0, convert(Float64, dist) / convert(Float64, instance.k), 1.0)
        blocks += ifelse(dist > 0, 0, 1)
        score += deceptive_dist
    end

    return score, blocks
end

"""
    bb_wrap_kdp(instance :: KDPInstance)

Wrap an instance of the Knjazew's Deceptive Permutation Problem and return a black box
that takes an assignment and returns its fitness value on the given instance.
"""
function bb_wrap_kdp(instance :: KDPInstance)
    # Preallocate memory
    lo :: Vector{Int64} = ones(Int64, instance.k)
    ip :: Vector{Int64} = zeros(Int64, instance.n)
    
    function evaluate(assignment :: Vector{Int64}) :: Float64
        return evaluate_kdp(instance, assignment, lo, ip)[1]
    end
    return evaluate
end

"""
    get_blocks_idx(insta :: KDPInstance)

Construct a list of blocks, groups of indices that are part of the same subfunction.
"""
function get_blocks_idx(insta :: KDPInstance) :: Vector{Vector{Int64}}
    idxs = zeros(Int64, insta.n)
    for (i, j) in enumerate(insta.mapping)
        idxs[j] = i
    end
    return [idxs[i*insta.k+1:(i+1)*insta.k] for i in 0:(insta.b-1)]
end