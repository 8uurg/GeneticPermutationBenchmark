# This file contains utilities for working with the Permutation Flowshop Problem
# such as Storage (Struct), Parsing, Evaluation and creating a Black-Box function. 
# !!!note The Permutation Flowshop Problem is a **minimization** problem.
#         As such the bb_wrap_pfs function negates the objective value
#         to allow maximization to work.

"A struct representing an instance of the Permutation Flowshop Problem"
struct PFSInstance{T <: Real}
    n :: Int64
    m :: Int64

    P :: Matrix{T} 
end

"""
    evaluate_pfs(instance :: PFSInstance{T}, permutation :: Vector{Int64};
        t_machine = zeros(T, instance.m))

Evaluate `permutation` against PFSInstance `instance`, and 
return its objective value: the flow time of the schedule.

The variable `t_machine` is the vector used to store the current completion time 
of each machine during the process of evaluating the schedule.

!!!note It is recommended to reuse the same vector `t_machine` for every evaluation.
        This avoids allocating and deallocating memory, and therefore relieves
        some pressure off the garbage collector.
"""
function evaluate_pfs(instance :: PFSInstance{T}, permutation :: Vector{Int64},
        t_machine :: Vector{T} = zeros(T, instance.m)) :: T where {T <: Int64}
    # Make sure the current completion time for each machine is zero.
    fill!(t_machine, zero(T))

    # Go over each job in the order of the permutation
    for j in permutation
        # The first task of the job always completes at start + P[1, j]
        t_machine[1] = t_machine[1] + instance.P[1, j]
        # The following tasks can only start when the machine is done 
        # /and/ the task at the previous machine is done.
        for x in 2:instance.m
            t_machine[x] = max(t_machine[x-1], t_machine[x]) + instance.P[x, j]
        end
    end

    # The final objective is the Makespan, or Flow Time.
    return maximum(t_machine)
end

"""
    bb_wrap_pfs(instance :: QAPInstance{T})

Wrap an instance of Permutation Flowshop Scheduling and return a black box
that takes a permutation and returns its fitness value on the given instance.
"""
function bb_wrap_pfs(instance :: PFSInstance{T}) where {T <: Real}
    # Preallocate the additional memory used in the evaluation.
    t_machine = zeros(T, instance.m)

    function evaluate(permutation :: Vector{Int64}) :: T
        return -evaluate_pfs(instance, permutation, t_machine)
    end

    return evaluate
end

"""
    parse_qap_qaplib(T, data :: String)

Parse a PFS instance from Taillard and return a QAPInstance struct. Argument `i` chooses
which instance to load.
"""
function parse_pfs_taillard(data :: String, i :: Int64) :: PFSInstance{Int64}
    # Select instance by index
    instance_separator = "number of jobs, number of machines, initial seed, upper bound and lower bound :"
    taillard_instance_strings = split(data, instance_separator, keepempty=false)
    # Only use the substring that intersts us.
    data = split(taillard_instance_strings[i], keepempty=false)
    n = parse(Int64, data[1])
    m = parse(Int64, data[2])
    # initial_seed = parse(Int64, data[3])
    # ub = parse(Int64, data[4])
    # lb = parse(Int64, data[5])
    # Elements 6-8 constitute ["processing", "times",":"]
    P = collect(reshape(parse.(Ref(Int64), data[9:9+n*m-1]), (n, m))')
    PFSInstance(n, m, P)
end
