struct NQSInstance
    n :: Int64
end

function evaluate_nqs(instance :: NQSInstance, solution :: Vector{Int64}, diag1_set :: BitVector, diag2_set :: BitVector)
    # Constants
    n = length(solution)
    @assert n == instance.n
    # Prepare memory structures
    resize!(diag1_set, n*2-1)
    resize!(diag2_set, n*2-1)
    fill!(diag1_set, false)
    fill!(diag2_set, false)
    #
    score = 0
    # 
    for (idx, e) in enumerate(solution)
        diag1_idx = idx - e + n
        diag2_idx = idx + e - 1
        if diag1_set[diag1_idx] || diag2_set[diag2_idx]
            score += 1
        end
        diag1_set[diag1_idx] = true
        diag2_set[diag2_idx] = true
    end
    return score
end

function bb_wrap_nqs(insta :: NQSInstance)
    diag1_set :: BitVector = falses(insta.n)
    diag2_set :: BitVector = falses(insta.n)
    function evaluate(perm :: Vector{Int64})
        return -evaluate_nqs(insta, perm, diag1_set, diag2_set)
    end
end