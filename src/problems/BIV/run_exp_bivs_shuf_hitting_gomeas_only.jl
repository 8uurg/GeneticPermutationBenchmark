println("Setting up...")
import Glob:glob
import DataFrames:DataFrame
using ProgressMeter
import CSV:CSV
using Random
# Note: Set JULIA_NUM_THREADS to the amount of threads to use.

# Number of runs, per approach, per instance
n_exp = 5
# (Maximum) amount of time for each run, per instance in seconds.
t_max = 100.0
# (Maximum) amount of evaluations
e_max = typemax(Int64) # 10000000
# Sidenote: An approach can converge and not use up the evaluations.

moments = [t_max]
moments_eval = [e_max]

if length(ARGS) > 0
    exp_idx_offset = parse(Int64, ARGS[1])
else
    exp_idx_offset = 0
end
path_results_time = "./results/results_exphit_biv_shuf_gomeas_$(exp_idx_offset)_time.csv"
path_results_evals = "./results/results_exphit_biv_shuf_gomeas_$(exp_idx_offset)_evals.csv"
path_results_time_hit = "./results/results_exphit_biv_shuf_gomeas_$(exp_idx_offset)_time_hit.csv"
path_results_evals_hit = "./results/results_exphit_biv_shuf_gomeas_$(exp_idx_offset)_evals_hit.csv"

# Make sure ./results exists
if !isdir("./results")
    mkdir("./results")
end

println("Loading problem & approaches...")
# Load problem utilities
include("./biv.jl")
# Load performance tracking utilities
include("../../utilities/trace.jl")
# Include permutation remapper
include("../permutationtools.jl")
# Load approaches
include("../../approaches/pGOMEA/pGOMEA.jl")
include("../../approaches/qGOMEA/qGOMEA.jl")
include("../../approaches/SimpleGA/IPSimpleGA.jl")
include("../../approaches/SimpleGA/RKSimpleGA.jl")

# And set them up
approaches = [
    # Permutation GOMEA
    # ("Permutation GOMEA - LT/Original - 1x FI", 
    #     (f, n, t, e; target_fitness) -> optimize_pgomea(f, n, t, e, forced_improvement = :original, target_fitness=target_fitness)),
    ("Permutation GOMEA - LT/Original - 10x FI", 
        (f, n, t, e; target_fitness) -> optimize_pgomea(f, n, t, e, forced_improvement = :extended, target_fitness=target_fitness)),
    # ("Permutation GOMEA - LT/Original - No FI", 
    #     (f, n, t, e; target_fitness) -> optimize_pgomea(f, n, t, e, forced_improvement = :none, target_fitness=target_fitness)),
    # ("Permutation GOMEA - RT - 1x FI", 
    #     (f, n, t, e; target_fitness) -> optimize_pgomea(f, n, t, e, forced_improvement = :original, fos_type=:random, target_fitness=target_fitness)),
    ("Permutation GOMEA - RT - 10x FI", 
        (f, n, t, e; target_fitness) -> optimize_pgomea(f, n, t, e, forced_improvement = :extended, fos_type=:random, target_fitness=target_fitness)),
    # ("Permutation GOMEA - RT - No FI", 
    #     (f, n, t, e; target_fitness) -> optimize_pgomea(f, n, t, e, forced_improvement = :none, fos_type=:random, target_fitness=target_fitness)),
    # ("Permutation GOMEA  - LT/Distance", 
    #     (f, n, t, e; target_fitness) -> optimize_pgomea(f, n, t, e, fos_type=:distance, target_fitness=target_fitness)),
    
    # qGOMEA
    # ("qGOMEA - LT/Distance - 1x FI - OX", 
    #     (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, forced_improvement = :original, target_fitness=target_fitness)),
    ("qGOMEA - LT/Distance - 10x FI - OX", 
        (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, forced_improvement = :extended, target_fitness=target_fitness)),
    ("qGOMEA - LT/Distance - No FI - OX", 
        (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, forced_improvement = :none, target_fitness=target_fitness)),
    
    # ("qGOMEA - LT/PermutationGOMEA Original - 10x FI - OX", 
    #     (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, fos_type=:original, target_fitness=target_fitness)),
    # ("qGOMEA - LT/PermutationGOMEA Original - No FI - OX", 
    #     (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, fos_type=:original, forced_improvement = :none, target_fitness=target_fitness)),
    # ("qGOMEA - LT/Distance - 10x FI - PMX", 
    #     (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, permutation_repair=:pmx, target_fitness=target_fitness)),

    # ("qGOMEA - RT - 1x FI - OX", 
    #     (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, forced_improvement = :original, fos_type=:random, target_fitness=target_fitness)),
    ("qGOMEA - RT - 10x FI - OX", 
        (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, forced_improvement = :extended, fos_type=:random, target_fitness=target_fitness)),
    ("qGOMEA - RT - No FI - OX", 
        (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, forced_improvement = :none, fos_type=:random, target_fitness=target_fitness)),
    
    # Random Key SimpleGA
    # ("Random Key SimpleGA", (f, n, t, e; target_fitness) -> optimize_rksimplega(f, n, t, e, target_fitness=target_fitness)),
    
    # Integer Permutation SimpleGA with various permutation crossover operators.
    # ("Integer Permutation SimpleGA - PMX", 
    #     (f, n, t, e; target_fitness) -> optimize_ipsimplega(PMX(n), f, n, t, e, target_fitness=target_fitness)),
    # ("Integer Permutation SimpleGA - OX", 
    #     (f, n, t, e; target_fitness) -> optimize_ipsimplega(OX(n), f, n, t, e, target_fitness=target_fitness)),
    # ("Integer Permutation SimpleGA - CX", 
    #     (f, n, t, e; target_fitness) -> optimize_ipsimplega(CX(n), f, n, t, e, target_fitness=target_fitness)),
    # ("Integer Permutation SimpleGA - ER", 
    #     (f, n, t, e; target_fitness) -> optimize_ipsimplega(ER(n), f, n, t, e, target_fitness=target_fitness)),
]

println("Initializing instances")
# Create instances.

instances = [
    # Inversions benchmark function
    # ("Inversion n=10" , BIVInstance(sorted_inversion,  10), convert(Float64, div( 10*( 10-1), 2))),
    # ("Inversion n=15" , BIVInstance(sorted_inversion,  15), convert(Float64, div( 15*( 15-1), 2))),
    # ("Inversion n=20" , BIVInstance(sorted_inversion,  20), convert(Float64, div( 20*( 20-1), 2))),
    # ("Inversion n=25" , BIVInstance(sorted_inversion,  25), convert(Float64, div( 25*( 25-1), 2))),
    # ("Inversion n=50" , BIVInstance(sorted_inversion,  50), convert(Float64, div( 50*( 50-1), 2))),
    # ("Inversion n=100", BIVInstance(sorted_inversion, 100), convert(Float64, div(100*(100-1), 2))),
    # ("Inversion n=200", BIVInstance(sorted_inversion, 200), convert(Float64, div(200*(200-1), 2))),
    # ("Inversion n=400", BIVInstance(sorted_inversion, 400), convert(Float64, div(400*(400-1), 2))),
    # ("Inversion n=800", BIVInstance(sorted_inversion, 800), convert(Float64, div(800*(800-1), 2))),
    # Sequential Inversions benchmark function
    # ("Sequential Inversion n=10" , BIVInstance(sorted_sequential_inversion,  10), convert(Float64,  10 - 1)),
    # ("Sequential Inversion n=15" , BIVInstance(sorted_sequential_inversion,  15), convert(Float64,  15 - 1)),
    # ("Sequential Inversion n=20" , BIVInstance(sorted_sequential_inversion,  20), convert(Float64,  20 - 1)),
    # ("Sequential Inversion n=25" , BIVInstance(sorted_sequential_inversion,  25), convert(Float64,  25 - 1)),
    # ("Sequential Inversion n=50" , BIVInstance(sorted_sequential_inversion,  50), convert(Float64,  50 - 1)),
    ("Sequential Inversion n=100", BIVInstance(sorted_sequential_inversion, 100), convert(Float64, 100 - 1)),
    ("Sequential Inversion n=200", BIVInstance(sorted_sequential_inversion, 200), convert(Float64, 200 - 1)),
    ("Sequential Inversion n=400", BIVInstance(sorted_sequential_inversion, 400), convert(Float64, 400 - 1)),
    ("Sequential Inversion n=800", BIVInstance(sorted_sequential_inversion, 800), convert(Float64, 800 - 1)),
    # Sequential Pairs benchmark function
    # ("Sequential Pairs n=10" , BIVInstance(sorted_sequential_pairs,  10), convert(Float64, 10 - 1)),
    # ("Sequential Pairs n=15" , BIVInstance(sorted_sequential_pairs,  15), convert(Float64, 15 - 1)),
    # ("Sequential Pairs n=20" , BIVInstance(sorted_sequential_pairs,  20), convert(Float64, 20 - 1)),
    # ("Sequential Pairs n=25" , BIVInstance(sorted_sequential_pairs,  25), convert(Float64, 25 - 1)),
    # ("Sequential Pairs n=50" , BIVInstance(sorted_sequential_pairs,  50), convert(Float64, 50 - 1)),
    ("Sequential Pairs n=100", BIVInstance(sorted_sequential_pairs, 100), convert(Float64,100 - 1)),
    ("Sequential Pairs n=200", BIVInstance(sorted_sequential_pairs, 200), convert(Float64,200 - 1)),
    ("Sequential Pairs n=400", BIVInstance(sorted_sequential_pairs, 400), convert(Float64,400 - 1)),
    ("Sequential Pairs n=800", BIVInstance(sorted_sequential_pairs, 800), convert(Float64,800 - 1)),
]

# Initialize results storage
results_time = DataFrame(
    [   String,    String,  Int64, Float64,    Float64], 
    [:instance, :approach, :exp_i,   :time, :objective])
results_eval = DataFrame(
    [   String,    String,  Int64,        Int64,    Float64], 
    [:instance, :approach, :exp_i, :evaluations, :objective])
results_time_hit = DataFrame(
    [   String,    String,  Int64, Float64,    Float64], 
    [:instance, :approach, :exp_i,   :time, :objective])
results_eval_hit = DataFrame(
    [   String,    String,  Int64,        Int64,    Float64], 
    [:instance, :approach, :exp_i, :evaluations, :objective])

begin
    println("Warming up approaches")
    instance_warmup = instances[1][2]
    bb_warmup = bb_wrap_biv(instance_warmup)
    for (approach_name, optimize_approach) in approaches
        # ~2.0s should be sufficient to get Julia to compile and run everything.
        # Alternatively, 5000 evaluations indicate that most things have been ran as well.
        optimize_approach(bb_warmup, instance_warmup.n, 2.0, 5000; target_fitness=nothing)
    end
    # results_lock = Threads.SpinLock()

    results_time_thread = [copy(results_time) for _ in 1:Threads.nthreads()]
    results_eval_thread = [copy(results_eval) for _ in 1:Threads.nthreads()]
    results_time_hit_thread = [copy(results_time_hit) for _ in 1:Threads.nthreads()]
    results_eval_hit_thread = [copy(results_eval_hit) for _ in 1:Threads.nthreads()]

    experiments = collect(Iterators.product(instances, approaches, 1:n_exp))
    
    println("Starting experiment, running a total of $(length(experiments)) experiments on $(Threads.nthreads()) thread(s).")

    progress = Progress(length(experiments))
    # @distributed
    Threads.@threads for ((instance_name, instance, instance_opt), (approach_name, optimize_approach), exp_i) in experiments
        
        bbf = random_remap(bb_wrap_biv(instance), instance.n)
        # Test evaluation.
        bbf(shuffle!(collect(1:instance.n)))
        # Perform GC for good measure.
        GC.gc()
        # Setup performance/time trace in black box.
        trace = ProgressTrace(moments, moments_eval, instance_opt)
        bbf_traced = trace_bb(bbf, trace)
        # Run experiment
        res = optimize_approach(bbf_traced, instance.n, t_max, e_max; target_fitness=instance_opt)
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
        push!(results_time_hit_thread[Threads.threadid()], (instance_name, approach_name, exp_i + exp_idx_offset, trace.hitting_time, instance_opt))
        push!(results_eval_hit_thread[Threads.threadid()], (instance_name, approach_name, exp_i + exp_idx_offset, trace.hitting_eval, instance_opt))
        next!(progress)
        # unlock(results_lock)
    end
    # Gather results for each thread.
    for th in 1:Threads.nthreads()
        append!(results_time, results_time_thread[th])
        append!(results_eval, results_eval_thread[th])
        append!(results_time_hit, results_time_hit_thread[th])
        append!(results_eval_hit, results_eval_hit_thread[th])
    end
    # Store results
    results_time |> CSV.write(path_results_time)
    results_eval |> CSV.write(path_results_evals)
    results_time_hit |> CSV.write(path_results_time_hit)
    results_eval_hit |> CSV.write(path_results_evals_hit)
end
