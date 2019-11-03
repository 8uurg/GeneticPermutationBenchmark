# This file contains utilities for working with the Quadratic Assignment Problem
# such as Storage (Struct), Parsing, Evaluation and creating a Black-Box function.
# !!!note The Quadratic Assignment Problem is a **minimization** problem.
#         As such the bb_wrap_qap function negates the objective value
#         to allow maximization to work.

"A struct representing an instance of the Quadratic Assignment Problem"
struct QAPInstance{T <: Real}
    n :: Int64

    A :: Matrix{T}
    B :: Matrix{T}
end

"""
    evaluate_qap(instance :: QAPInstance{T}, assignment :: Vector{Int64})

Evaluate `assignment` against QAPInstance `instance`, and 
return its objective value.
"""
function evaluate_qap(instance :: QAPInstance{T}, 
        assignment :: Vector{Int64}) :: T where {T <: Real}
    # Initially score is zero.
    score :: T = zero(T)

    # ∑i∑j Aᵢⱼ ⋅ Bπᵢπⱼ
    @inbounds for i in 1:instance.n, j in 1:instance.n 
        score += instance.A[i, j] * instance.B[assignment[i], assignment[j]] 
    end

    return score
end

"""
    bb_wrap_qap(instance :: QAPInstance{T})

Wrap an instance of the Quadradic Assignment Problem and return a black box
that takes an assignment and returns its fitness value on the given instance.
"""
function bb_wrap_qap(instance :: QAPInstance{T}) where {T <: Real}
    function evaluate(assignment :: Vector{Int64}) :: T
        return -evaluate_qap(instance, assignment)
    end
    return evaluate
end

"""
    parse_qap_qaplib(T, data :: String)

Parse a QAP instance from QAPLIB and return a QAPInstance struct.
"""
function parse_qap_qaplib(T, data :: String) :: QAPInstance{T}
    # Split on newlines, spaces, do not keep empty strings.
    data = split(data, keepempty=false)
    # First item is the amount of 'facilities' and 'locations'.
    n = parse(Int64, data[1])
    # The matrices are consitutent of the remaining elements
    # containing n*n each. 
    A = collect(reshape(parse.(Ref(T), data[2:1+n*n]), (n, n))')
    B = collect(reshape(parse.(Ref(T), data[2+n*n:1+2*n*n]), (n, n))')
    QAPInstance{T}(n, A, B)
end

"""
Relabel the A matrix of a Quadratic Assignment Problem instance.

Note: this changes the meaning of a solution. 
A solution `s` in for the unrelabeled `instance` is permutated by the the assignment's.
hence `s[assignment]` yields a solution with the same objective value.
"""
function relabel_A(instance :: QAPInstance, assignment :: Vector{Int64})
    QAPInstance(instance.n, instance.A[assignment, assignment], instance.B)
end

"""
Relabel the B matrix of a Quadratic Assignment Problem instance.

Note: this changes the meaning of a solution.
A solution `s` in for the unrelabeled `instance` is permutated by the the assignment's inverse.
hence `s[invperm(assignment)]` yields a solution with the same objective value.
If `assignment` was the optimal solution then the new optimal solution is `[1, ..., instance.n]`.
"""
function relabel_B(instance :: QAPInstance, assignment :: Vector{Int64})
    QAPInstance(instance.n, instance.A, instance.B[assignment, assignment])
end
 
"""
Relabel the A matrix of a Quadratic Assignment Problem instance.

Note: this changes the meaning of a solution. 
A solution `s` in for the unrelabeled `instance` is permutated by the the assignment's inverse.
If `assignment` was the optimal solution then
the new optimal solution is `[1, ..., instance.n]`.
"""
function relabel_inv_A(instance :: QAPInstance, assignment :: Vector{Int64})
    iassignment = invperm(assignment)
    QAPInstance(instance.n, instance.A[iassignment, iassignment], copy(instance.B))
end

"""
Relabel the B matrix of a Quadratic Assignment Problem instance.

Note: this changes the meaning of a solution.
A solution `s` in for the unrelabeled `instance` is permutated by the the assignment's.

"""
function relabel_inv_B(instance :: QAPInstance, assignment :: Vector{Int64})
    iassignment = invperm(assignment)
    QAPInstance(instance.n, copy(instance.A), instance.B[iassignment, iassignment])
end

"""
    qap_distance_weight_matrix(n :: Int64)

Create a symmetric n × n matrix where each position equals the distance
between the two positions. 
"""
function qap_distance_weight_matrix(n :: Int64)
    m = zeros(Int64, (n, n))
    for r in 1:n, c in 1:n
        m[r, c] = abs(r - c)
    end
    return m
end


"""
Create a symmetric n × n matrix where each position equals n - the distance
between the two positions. 
"""
function qap_inverted_distance_weight_matrix(n :: Int64)
    m = zeros(Int64, (n, n))
    for r in 1:n, c in 1:n
        m[r, c] = n - abs(r - c)
    end
    return m
end

"""
    qap_bidirectional_path_weight_matrix(n :: Int64)

Construct a weight matrix corresponding to finding the shortest (bidirectional)
Hamiltonian path (a path that visits all vertices exactly once).
"""
function qap_bidirectional_path_weight_matrix(n :: Int64)
    m = zeros(Int64, (n, n))
    for r in 1:(n-1)
        m[r, r+1] = 1
        m[r+1, r] = 1
    end
    return m
end

"""
    qap_tsp_weight_matrix(n :: Int64)

Construct a weight matrix corresponding to the Travelling Salesperson Problem
"""
function qap_tsp_weight_matrix(n :: Int64)
    m = zeros(Int64, (n, n))
    for r in 1:(n-1)
        m[r, r+1] = 1
    end
    m[n, 1] = 1
    return m
end

# p = optimize_gomea(bb_wrap_qap(QAPInstance(insta.n, qap_inverted_distance_weight_matrix(insta.n), copy(insta.A))), insta.n, 10.0; forced_improvement=:extended)