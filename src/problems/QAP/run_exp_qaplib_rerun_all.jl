println("Setting up...")
import Glob:glob
import DataFrames:DataFrame
using ProgressMeter
import CSV:CSV
using Random

# Number of runs, per approach, per instance
n_exp = 12
# (Maximum) amount of time for each run, per instance in seconds.
# For this experiment: 10 minutes.
t_max = 10.0 * 60.0
# (Maximum) amount of evaluations
e_max = typemax(Int64) # 10000000
# Sidenote: An approach can converge and not use up the evaluations.

moments = [t_max]
moments_eval = [e_max]

# if length(ARGS) > 0
#     exp_idx_offset = parse(Int64, ARGS[1])
# else
approach_idx = parse(Int64, ARGS[1])
exp_idx_offset = rand(UInt32)
# end
path_results_time = "./results/results_qap_qaplib_rrg_$(approach_idx)_$(exp_idx_offset)_time.csv"
path_results_evals = "./results/results_qap_qaplib_rrg_$(approach_idx)_$(exp_idx_offset)_evals.csv"
path_results_hittime = "./results/results_qap_qaplib_rrg_$(approach_idx)_$(exp_idx_offset)_hittime.csv"
path_results_hitevals = "./results/results_qap_qaplib_rrg_$(approach_idx)_$(exp_idx_offset)_hitevals.csv"
path_instances = "./instances/qaplib/instances"

# Make sure ./results exists
if !isdir("./results")
    mkdir("./results")
end

# Find instance files.
instance_names = Set(string.(["chr22b", "els19", "esc32b", "esc128", "kra30b", "lipa30b", 
    "lipa90a", "lipa90b", "nug30", "rou20", "scr15", "scr20", "sko64", "ste36c", 
    "tai35a", "tai100b", "tai256c", "tho40", "tho150", "wil50"], ".dat"))
# instance_names = Set(string.(["chr22b", "els19", "esc32b"], ".dat"))
instance_bounds = CSV.read("./instances/qaplib/bounds.txt")
instance_bounds[!, :instance] = string.(lowercase.(instance_bounds[!, :name]), ".dat")
instances = filter(row -> (row[:instance] in instance_names), instance_bounds)
instances[!, :path] = string.(path_instances, "/", instances[!, :instance])

# Check if instances are present.
err_no_instances = """Not all 20 instances found."""
@assert size(instances, 1) == length(instance_names) err_no_instances

println("Loading problem & approaches...")
# Load problem utilities
include("./qap.jl")
# Load performance tracking utilities
include("../../utilities/trace.jl")
# Load approaches
# Note: Do not load the original PermutationGOMEA code, as it requires a compiled binary.
# include("../../approaches/PermutationGOMEA/PermutationGOMEA.jl")
include("../../approaches/pGOMEA/pGOMEA.jl")
include("../../approaches/qGOMEA/qGOMEA.jl")
include("../../approaches/SimpleGA/IPSimpleGA.jl")
include("../../approaches/SimpleGA/RKSimpleGA.jl")

# And set them up
approaches_under_test = [
    # Permutation GOMEA
    # ("Permutation GOMEA - LT/Original - 1x FI", 
    #     (f, n, t, e; target_fitness) -> optimize_pgomea(f, n, t, e; target_fitness=target_fitness, forced_improvement = :original)),
    ("Permutation GOMEA - LT/Original - 10x FI", 
        (f, n, t, e; target_fitness) -> optimize_pgomea(f, n, t, e; target_fitness=target_fitness, forced_improvement = :extended)),
    # ("Permutation GOMEA - LT/Original - No FI", 
    #     (f, n, t, e; target_fitness) -> optimize_pgomea(f, n, t, e; target_fitness=target_fitness, forced_improvement = :none)),
    # ("Permutation GOMEA - RT - 1x FI", 
    #     (f, n, t, e; target_fitness) -> optimize_pgomea(f, n, t, e; target_fitness=target_fitness, forced_improvement = :original, fos_type=:random)),
    ("Permutation GOMEA - RT - 10x FI", 
        (f, n, t, e; target_fitness) -> optimize_pgomea(f, n, t, e; target_fitness=target_fitness, forced_improvement = :extended, fos_type=:random)),
    # ("Permutation GOMEA - RT - No FI", 
    #     (f, n, t, e; target_fitness) -> optimize_pgomea(f, n, t, e; target_fitness=target_fitness, forced_improvement = :none, fos_type=:random)),
    # ("Permutation GOMEA  - LT/Distance", 
    #     (f, n, t, e; target_fitness) -> optimize_pgomea(f, n, t, e; target_fitness=target_fitness, fos_type=:distance)),
    
    # qGOMEA
    # ("qGOMEA - LT/Distance - 1x FI - OX", 
    #     (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e; target_fitness=target_fitness, forced_improvement = :original)),
    # ("qGOMEA - LT/Distance - 10x FI - OX", 
    #     (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e; target_fitness=target_fitness, forced_improvement = :extended)),
    ("qGOMEA - LT/Distance - No FI - OX", 
        (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e; target_fitness=target_fitness, forced_improvement = :none)),
    
    # ("qGOMEA - LT/PermutationGOMEA Original", 
    #     (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e; target_fitness=target_fitness, fos_type=:original)),
    # ("qGOMEA - LT/Distance - 10x FI - PMX", 
    #     (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e; target_fitness=target_fitness, permutation_repair=:pmx)),
    ("qGOMEA - LT/Distance - No FI - PMX", 
        (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e; target_fitness=target_fitness, forced_improvement = :none, permutation_repair=:pmx)),
    ("qGOMEA - RT - No FI - PMX", 
        (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e; target_fitness=target_fitness, forced_improvement = :none, permutation_repair=:pmx, fos_type=:random)),


    # ("qGOMEA - RT - 1x FI - OX", 
    #     (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e; target_fitness=target_fitness, forced_improvement = :original, fos_type=:random)),
    # ("qGOMEA - RT - 10x FI - OX", 
    #     (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e; target_fitness=target_fitness, forced_improvement = :extended, fos_type=:random)),
    ("qGOMEA - RT - No FI - OX", 
        (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e; target_fitness=target_fitness, forced_improvement = :none, fos_type=:random)),
    
    # Random Key SimpleGA
    ("Random Key SimpleGA", (f, n, t, e; target_fitness) -> optimize_rksimplega(f, n, t, e; target_fitness=target_fitness)),
    
    # Integer Permutation SimpleGA with various permutation crossover operators.
    ("Integer Permutation SimpleGA - PMX", 
        (f, n, t, e; target_fitness) -> optimize_ipsimplega(PMX(n), f, n, t, e; target_fitness=target_fitness)),
    ("Integer Permutation SimpleGA - OX", 
        (f, n, t, e; target_fitness) -> optimize_ipsimplega(OX(n), f, n, t, e; target_fitness=target_fitness)),
    ("Integer Permutation SimpleGA - CX", 
        (f, n, t, e; target_fitness) -> optimize_ipsimplega(CX(n), f, n, t, e; target_fitness=target_fitness)),
    ("Integer Permutation SimpleGA - ER", 
        (f, n, t, e; target_fitness) -> optimize_ipsimplega(ER(n), f, n, t, e; target_fitness=target_fitness)),

    # Original PermutationGOMEA
    # ("Original Implementation Permutation GOMEA", (f, n, t, e; target_fitness) -> optimize_permutationgomea(f, n, t, e)),
]

if approach_idx <= 0
    error("Index of approach under test is zero or negative: invalid value.")
end

if approach_idx > length(approaches_under_test)
    error("Index of approach under test is greater than the number available: invalid value.")
end

approaches = [
    approaches_under_test[approach_idx]
]

println("Reading and parsing instances")
# Load and parse all instances.
instances = collect(Iterators.flatten(
    begin 
        instance_file = open(instance_row[:path], "r"); 
        instance_str = read(instance_file, String);
        close(instance_file);
        ( (instance_row[:instance], instance, instance_row[:lb], instance_row[:ub]) 
            for (i, instance) in 
                enumerate([parse_qap_qaplib(Int64, instance_str)]))
    end 
    for instance_row in eachrow(instances)))

# Initialize results storage
results_time = DataFrame(
    [Vector{String}(), Vector{String}(), Vector{Int64}(), Vector{Float64}(), Vector{Float64}()], 
    [       :instance,        :approach,          :exp_i,             :time,        :objective])
results_eval = DataFrame(
    [Vector{String}(), Vector{String}(), Vector{Int64}(), Vector{Int64}(), Vector{Float64}()], 
    [       :instance,        :approach,          :exp_i,    :evaluations,        :objective])

results_hittime = DataFrame(
    [Vector{String}(), Vector{String}(), Vector{Int64}(), Vector{Float64}(),    Vector{Float64}()], 
    [       :instance,        :approach,          :exp_i,             :time,           :objective])
results_hiteval = DataFrame(
    [Vector{String}(), Vector{String}(), Vector{Int64}(), Vector{Int64}(), Vector{Float64}()], 
    [       :instance,        :approach,          :exp_i,    :evaluations,        :objective])

begin
    println("Warming up approaches")
    instance_warmup = instances[1][2]
    bb_warmup = bb_wrap_qap(instance_warmup)
    for (approach_name, optimize_approach) in approaches
        # ~2.0s should be sufficient to get Julia to compile and run everything.
        # Alternatively, 5000 evaluations indicate that most things have been ran as well.
        trace = ProgressTrace(moments, moments_eval, 0.0)
        bb_warmup_traced = trace_bb(bb_warmup, trace)
        optimize_approach(bb_warmup_traced, instance_warmup.n, 2.0, 5000; target_fitness=nothing)
        postprocess_trace!(trace)
    end
    # results_lock = Threads.SpinLock()

    experiments = collect(Iterators.product(instances, approaches, 1:n_exp))
    
    println("Starting experiment, running a total of $(length(experiments)) experiments.")

    progress = Progress(length(experiments))
    # @distributed
    for ((instance_name, instance, lb, ub), (approach_name, optimize_approach), exp_i) in experiments
        
        bbf = bb_wrap_qap(instance)
        # Test evaluation.
        bbf(shuffle!(collect(1:instance.n)))
        # Perform GC for good measure.
        GC.gc()
        lb_target = convert(Float64, -lb)
        # Setup performance/time trace in black box.
        trace = ProgressTrace(moments, moments_eval, lb_target)
        bbf_traced = trace_bb(bbf, trace)
        # Run experiment
        res = optimize_approach(bbf_traced, instance.n, t_max, e_max, target_fitness=lb_target)
        # Postprocess trace (eg. clean it up)
        postprocess_trace!(trace)
        # Dump results.
        # lock(results_lock)
        for (time, fitness) in zip(trace.moments, trace.results)
            push!(results_time, (instance_name, approach_name, exp_i + exp_idx_offset, time, -fitness))
        end
        for (evals, fitness) in zip(trace.moments_n_eval, trace.results_eval)
            push!(results_eval, (instance_name, approach_name, exp_i + exp_idx_offset, evals, -fitness))
        end
        push!(results_hittime, (instance_name, approach_name, exp_i + exp_idx_offset, trace.hitting_time, -trace.best))
        push!(results_hiteval, (instance_name, approach_name, exp_i + exp_idx_offset, trace.hitting_eval, -trace.best))

        next!(progress)
        # unlock(results_lock)
    end
    # Store results
    results_time |> CSV.write(path_results_time)
    results_eval |> CSV.write(path_results_evals)
    results_hittime |> CSV.write(path_results_hittime)
    results_hiteval |> CSV.write(path_results_hitevals)
end
