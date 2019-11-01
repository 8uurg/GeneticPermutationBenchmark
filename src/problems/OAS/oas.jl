# This file contains utilities for working with the Order Acceptance and Scheduling with Sequence Dependent Setup Times Problem
# such as Storage (Struct), Parsing, Evaluation and creating a Black-Box function. 
import OffsetArrays:OffsetArray

# Definition for Order Acceptance and Scheduling (hereafter OAS)
# Written according to the definition in
# Cesaret, Bahriye, Ceyda Oǧuz, and F Sibel Salman (2012).
#  “A tabu search algorithm for order acceptance and scheduling”.
#  In: Computers & Operations Re- search 39.6, pp. 1197–1205.

"A struct representing an instance of the Order Acceptance and Scheduling with Sequence Dependent Setup Times Problem"
struct OASInstance
    # Amount of orders
    n :: Int64

    # An instance is defined by, for each order i
    # The release times r_i: the time at which an order can be performed
    r :: Array{Int64, 1}
    # The processing time p_i: the time it takes to perform an order
    p :: Array{Int64, 1}
    # The due date d_i: the point after which a penalty is incurred for every point in time.
    d :: Array{Int64, 1}
    # The deadline d̄_i: the point in time after which no revenue can be gained anymore.
    d̄ :: Array{Int64, 1}
    # The maximum revenue w_i
    e :: Array{Int64, 1}
    # The weight (tardiness penalty after the due date) w_i
    # Is actually set to e ./ (d̄ .- d), as this ensures that the deadline is the point
    # in time in which no more revenue can be gained.
    w :: Array{Float64, 1}
    # For every order i j, a sequence dependent setup time s_ij
    # Also, a special case i=0 for no preceding order s_0j
    # Uses OffsetArrays to make indexing with 0 possible.
    s :: OffsetArray{Int64, 2, Array{Int64, 2}}
end

"""
    evaluate_oas(instance :: OASInstance, solution :: Vector{Int64})
    evaluate_oas(instance :: OASInstance, solution :: Vector{Int64}; selected :: BitVector)

A fast implementation for evaluating a permutation `solution` against an OASInstance `instance`,
returning the profit obtained, as well as the completion time of the last order.

Skips over orders that are past their deadline, rather than returning a value corresponding to an infeasible solution.

The variant with selected as named argument will skip over items with their bit turned off,
as well as turn off the bits of unselected orders.
"""
function evaluate_oas(instance :: OASInstance, solution :: Vector{Int64}; selected :: Union{Nothing, BitVector}=nothing)
    n = length(solution)
    profit :: Float64 = 0.0
    Ctemp :: Int64 = typemin(Int64)
    solution′ᵢ = 0
    for i in 1:n
        @inbounds solutionᵢ = solution[i]
        # Note solution′ᵢ can be zero, make sure instance.s can be indexed for the first
        # axis with 0 by using - for example - an OffsetArray

        # Skip unselected items.
        if selected != nothing && !selected[solutionᵢ]
            continue
        end

        Ctemp′ = Ctemp
        @inbounds Ctemp = max(Ctemp, instance.r[solutionᵢ]) + instance.s[solution′ᵢ, solutionᵢ] + instance.p[solutionᵢ]
        # This is in their algorithm, but I am pretty sure this is wrong.
        # if Ctemp > instance.d̄[solution[a]]
        # Rather, the text seems to imply breaking at the first item
        # That violates the deadline.
        @inbounds if Ctemp >= instance.d̄[solutionᵢ]
            Ctemp = Ctemp′
            if selected != nothing
                selected[solutionᵢ] = false
            end
            #solution′ᵢ = solutionᵢ
            continue
        end
        @inbounds if Ctemp > instance.d[solutionᵢ]
            @inbounds profit += instance.e[solutionᵢ] + instance.w[solutionᵢ] * (instance.d[solutionᵢ] - Ctemp)
        else
            @inbounds profit += instance.e[solutionᵢ]
        end
        solution′ᵢ = solutionᵢ
    end
    # Return typemax for pastd̄i, so that we can take the minimum over blocks
    # to indicate the first index that goes past a deadline.
    return profit, Ctemp
end

"""
bb_wrap_oas(instance :: OASInstance)

Wrap an instance of the Order Acceptance and Scheduling [...] Problem and return a black box
that takes an permutation and returns its fitness value on the given instance.
"""
function bb_wrap_oas(instance :: OASInstance)

    function evaluate(permutation :: Vector{Int64})
        profit, Ctemp = evaluate_oas(instance, permutation)
        return profit
    end

    return evaluate
end

"""
    invperm!(dst :: Vector{Int64}, src :: Vector{Int64})

Write the inverse permutation to dst -- rather than storing it in a newly allocated piece of memory.

!!!warning This function does not perform bound-checks, nor does it resize the destination vector. 
"""
function invperm!(dst :: Vector{Int64}, src :: Vector{Int64})
    @inbounds for (i, j) in enumerate(src)
        dst[j] = i
    end
    dst
end


"""
bb_wrap_oas_assignment(instance :: OASInstance)

Wrap an instance of the Order Acceptance and Scheduling [...] Problem and return a black box
that takes an assignment and returns its fitness value on the given instance.
"""
function bb_wrap_oas_assignment(instance :: OASInstance)
    permutation = collect(1:instance.n)

    function evaluate(assignment :: Vector{Int64})
        invperm!(permutation, assignment)
        profit, Ctemp = evaluate_oas(instance, permutation)
        return profit
    end

    return evaluate
end
