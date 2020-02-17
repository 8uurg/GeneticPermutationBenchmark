println("Setting up...")
import Glob:glob
import DataFrames:DataFrame
using ProgressMeter
import CSV:CSV
using Random
# Note: Set JULIA_NUM_THREADS to the amount of threads to use.

# Number of runs, per approach, per instance
n_exp = 3
# (Maximum) amount of time for each run, per instance in seconds.
# (Maximum) amount of evaluations
e_max = typemax(Int64) # 10000000
# Sidenote: An approach can converge and not use up the evaluations.

moments_eval = [e_max]

if length(ARGS) > 0
    exp_idx_offset = parse(Int64, ARGS[1])
else
    exp_idx_offset = 0
end
path_results_time = "./results/results_oas_tm_p1_ces_$(exp_idx_offset)_time.csv"
path_results_evals = "./results/results_oas_tm_p1_ces_$(exp_idx_offset)_evals.csv"
path_instances = "./instances/Dataset_OAS"

# Make sure ./results exists
if !isdir("./results")
    mkdir("./results")
end

# Find instance files.
instances = glob(string(path_instances, "/*/*/*/*.txt"))

# Check if instances are present.
err_no_instances = """No instances found."""
@assert length(instances) != 0 err_no_instances

println("Loading problem & approaches...")
# Load problem utilities
include("./oas.jl")
# Load performance tracking utilities
include("../../utilities/trace.jl")
# Load approaches
include("../../approaches/pGOMEA/pGOMEA.jl")
include("../../approaches/qGOMEA/qGOMEA.jl")
include("../../approaches/SimpleGA/IPSimpleGA.jl")
include("../../approaches/SimpleGA/RKSimpleGA.jl")

# And set them up
approaches = [
    # Permutation GOMEA
    # ("Permutation GOMEA - LT/Original - 1x FI", 
    #     (f, n, t, e) -> optimize_pgomea(f, n, t, e, forced_improvement = :original)),
    ("Permutation GOMEA - LT/Original - 10x FI", 
        (f, n, t, e) -> optimize_pgomea(f, n, t, e, forced_improvement = :extended)),
    # ("Permutation GOMEA - LT/Original - No FI", 
    #     (f, n, t, e) -> optimize_pgomea(f, n, t, e, forced_improvement = :none)),
    # ("Permutation GOMEA - RT - 1x FI", 
    #     (f, n, t, e) -> optimize_pgomea(f, n, t, e, forced_improvement = :original, fos_type=:random)),
    ("Permutation GOMEA - RT - 10x FI", 
        (f, n, t, e) -> optimize_pgomea(f, n, t, e, forced_improvement = :extended, fos_type=:random)),
    # ("Permutation GOMEA - RT - No FI", 
    #     (f, n, t, e) -> optimize_pgomea(f, n, t, e, forced_improvement = :none, fos_type=:random)),
    # ("Permutation GOMEA  - LT/Distance", 
    #     (f, n, t, e) -> optimize_pgomea(f, n, t, e, fos_type=:distance)),
    
    # qGOMEA
    # ("qGOMEA - LT/Distance - 1x FI - OX", 
    #     (f, n, t, e) -> optimize_qgomea(f, n, t, e, forced_improvement = :original)),
    ("qGOMEA - LT/Distance - 10x FI - OX", 
        (f, n, t, e) -> optimize_qgomea(f, n, t, e, forced_improvement = :extended)),
    # ("qGOMEA - LT/Distance - No FI - OX", 
    #     (f, n, t, e) -> optimize_qgomea(f, n, t, e, forced_improvement = :none)),
    
    # ("qGOMEA - LT/PermutationGOMEA Original", 
    #     (f, n, t, e) -> optimize_qgomea(f, n, t, e, fos_type=:original)),
    # ("qGOMEA - LT/Distance - 10x FI - PMX", 
    #     (f, n, t, e) -> optimize_qgomea(f, n, t, e, permutation_repair=:pmx)),

    # ("qGOMEA - RT - 1x FI - OX", 
    #     (f, n, t, e) -> optimize_qgomea(f, n, t, e, forced_improvement = :original, fos_type=:random)),
    ("qGOMEA - RT - 10x FI - OX", 
        (f, n, t, e) -> optimize_qgomea(f, n, t, e, forced_improvement = :extended, fos_type=:random)),
    # ("qGOMEA - RT - No FI - OX", 
    #     (f, n, t, e) -> optimize_qgomea(f, n, t, e, forced_improvement = :none, fos_type=:random)),
    
    # Random Key SimpleGA
    # ("Random Key SimpleGA", (f, n, t, e) -> optimize_rksimplega(f, n, t, e)),
    
    # Integer Permutation SimpleGA with various permutation crossover operators.
    # ("Integer Permutation SimpleGA - PMX", 
    #     (f, n, t, e) -> optimize_ipsimplega(PMX(n), f, n, t, e)),
    # ("Integer Permutation SimpleGA - OX", 
    #     (f, n, t, e) -> optimize_ipsimplega(OX(n), f, n, t, e)),
    # ("Integer Permutation SimpleGA - CX", 
    #     (f, n, t, e) -> optimize_ipsimplega(CX(n), f, n, t, e)),
    # ("Integer Permutation SimpleGA - ER", 
    #     (f, n, t, e) -> optimize_ipsimplega(ER(n), f, n, t, e)),
]

println("Reading and parsing instances")
# Load and parse all instances.
instances = collect(Iterators.flatten(
    begin
        ( (basename(instance_path), instance) 
            for (i, instance) in 
                enumerate([parse_oas_cesaret(instance_path)]))
    end 
    for instance_path in instances))

# Initialize results storage
results_time = DataFrame(
    [   String,    String,  Int64, Float64,    Float64], 
    [:instance, :approach, :exp_i,   :time, :objective])
results_eval = DataFrame(
    [   String,    String,  Int64,        Int64,    Float64], 
    [:instance, :approach, :exp_i, :evaluations, :objective])

begin
    println("Warming up approaches")
    instance_warmup = instances[1][2]
    bb_warmup = bb_wrap_oas(instance_warmup)
    for (approach_name, optimize_approach) in approaches
        # ~2.0s should be sufficient to get Julia to compile and run everything.
        # Alternatively, 5000 evaluations indicate that most things have been ran as well.
        optimize_approach(bb_warmup, instance_warmup.n, 2.0, 5000)
    end
    # results_lock = Threads.SpinLock()

    results_time_thread = [copy(results_time) for _ in 1:Threads.nthreads()]
    results_eval_thread = [copy(results_eval) for _ in 1:Threads.nthreads()]

    shuffle!(instances)
    experiments = collect(Iterators.product(instances, approaches, 1:n_exp))

    println("Starting experiment, running a total of $(length(experiments)) experiments on $(Threads.nthreads()) thread(s).")

    progress = Progress(length(experiments))
    # @distributed
    Threads.@threads for ((instance_name, instance), (approach_name, optimize_approach), exp_i) in experiments
        t_max = convert(Float64, instance.n^2)/100+0.1
        moments = [t_max]
        bbf = bb_wrap_oas(instance)
        # Test evaluation.
        bbf(shuffle!(collect(1:instance.n)))
        # Perform GC for good measure.
        # GC.gc()
        # Setup performance/time trace in black box.
        trace = ProgressTrace(moments, moments_eval, 0.0)
        bbf_traced = trace_bb(bbf, trace)
        # Run experiment
        res = optimize_approach(bbf_traced, instance.n, t_max, e_max)
        # Postprocess trace (eg. clean it up)
        postprocess_trace!(trace)
        # Dump results.
        # lock(results_lock)
        for (time, fitness) in zip(trace.moments, trace.results)
            push!(results_time_thread[Threads.threadid()], (instance_name, approach_name, exp_i + exp_idx_offset, time, fitness))
        end
        for (evals, fitness) in zip(trace.moments_n_eval, trace.results_eval)
            push!(results_eval_thread[Threads.threadid()], (instance_name, approach_name, exp_i + exp_idx_offset, evals, fitness))
        end
        next!(progress)
        # unlock(results_lock)
    end
    # Gather results for each thread.
    for th in 1:Threads.nthreads()
        append!(results_time, results_time_thread[th])
        append!(results_eval, results_eval_thread[th])
    end
    # Store results
    results_time |> CSV.write(path_results_time)
    results_eval |> CSV.write(path_results_evals)
end
