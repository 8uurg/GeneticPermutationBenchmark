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
