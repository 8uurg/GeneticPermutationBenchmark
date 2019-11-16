println("Setting up...")
import Glob:glob
import DataFrames:DataFrame
using ProgressMeter
import CSV:CSV
using Random
# Note: Set JULIA_NUM_THREADS to the amount of threads to use.

# Number of runs, per approach, per instance
n_exp = 1#10
# (Maximum) amount of time for each run, per instance.
t_max = 1.0# Inf64 # 10.0
# (Maximum) amount of evaluations
e_max = typemax(Int64) # 10000000

moments = [t_max]
moments_eval = [e_max]

if length(ARGS) > 0
    exp_idx_offset = parse(Int64, ARGS[1])
else
    exp_idx_offset = 0
end
path_results_time = "./results/results_$(exp_idx_offset)_time.csv"
path_results_evals = "./results/results_$(exp_idx_offset)_evals.csv"
path_instances = "./instances/taillard/instances"

# Make sure ./results exists
if !isdir("./results")
    mkdir("./results")
end

# Find instance files.
instances = glob(string(path_instances, "/*.txt"))

# Check if instances are present.
err_no_instances = """No instances found."""
@assert length(instances) != 0 err_no_instances

println("Loading problem & approaches...")
# Load problem utilities
include("./pfs.jl")
# Load performance tracking utilities
include("../../utilities/trace.jl")
# Load approaches
include("../../approaches/pGOMEA/pGOMEA.jl")
include("../../approaches/qGOMEA/qGOMEA.jl")
include("../../approaches/SimpleGA/IPSimpleGA.jl")
include("../../approaches/SimpleGA/RKSimpleGA.jl")

# And set them up
approaches = [
    ("Permutation GOMEA", (f, n, t, e) -> optimize_pgomea(f, n, t, e)),
    ("qGOMEA", (f, n, t, e) -> optimize_qgomea(f, n, t, e)),
    ("Random Key SimpleGA", (f, n, t, e) -> optimize_rksimplega(f, n, t, e)),
    ("Integer Permutation (PMX) SimpleGA", 
        (f, n, t, e) -> optimize_ipsimplega(PMX(n), f, n, t, e)),
    ("Integer Permutation (OX) SimpleGA", 
        (f, n, t, e) -> optimize_ipsimplega(OX(n), f, n, t, e)),
    ("Integer Permutation (CX) SimpleGA", 
        (f, n, t, e) -> optimize_ipsimplega(CX(n), f, n, t, e)),
    ("Integer Permutation (ER) SimpleGA", 
        (f, n, t, e) -> optimize_ipsimplega(ER(n), f, n, t, e)),
]

println("Reading and parsing instances")
# Load and parse all instances.
instances = collect(Iterators.flatten(
    begin 
        instance_file = open(instance_path, "r"); 
        instance_str = read(instance_file, String);
        close(instance_file);
        ( (basename(instance_path), instance) 
            for (i, instance) in 
                enumerate(parse_pfs_taillard_all(instance_str)))
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
    bb_warmup = bb_wrap_pfs(instance_warmup)
    for (approach_name, optimize_approach) in approaches
        # ~2.0s should be sufficient to get Julia to compile and run everything.
        # Alternatively, 5000 evaluations indicate that most things have been ran as well.
        optimize_approach(bb_warmup, instance_warmup.n, 2.0, 5000)
    end
    results_lock = Threads.SpinLock()

    experiments = collect(Iterators.product(instances, approaches, 1:n_exp))
    
    println("Starting experiment, running a total of $(length(experiments)) experiments on $(Threads.nthreads()) thread(s).")

    progress = Progress(length(experiments))
    # @distributed
    Threads.@threads for ((instance_name, instance), (approach_name, optimize_approach), exp_i) in experiments
        
        bbf = bb_wrap_pfs(instance)
        # Test evaluation.
        bbf(shuffle!(collect(1:instance.n)))
        # Perform GC for good measure.
        GC.gc()
        # Setup performance/time trace in black box.
        trace = ProgressTrace(moments, moments_eval, 0.0)
        bbf_traced = trace_bb(bbf, trace)
        # Run experiment
        res = optimize_approach(bbf_traced, instance.n, t_max, e_max)
        # Postprocess trace (eg. clean it up)
        postprocess_trace!(trace)
        # Dump results.
        lock(results_lock)
        for (time, fitness) in zip(trace.moments, trace.results)
            push!(results_time, (instance_name, approach_name, exp_i + exp_idx_offset, time, -fitness))
        end
        for (evals, fitness) in zip(trace.moments_n_eval, trace.results_eval)
            push!(results_eval, (instance_name, approach_name, exp_i + exp_idx_offset, evals, -fitness))
        end
        next!(progress)
        unlock(results_lock)
    end
    # Store results
    results_time |> CSV.write(path_results_time)
    results_eval |> CSV.write(path_results_evals)
end