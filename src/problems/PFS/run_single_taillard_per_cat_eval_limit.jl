println("Setting up...")
import Glob:glob
import DataFrames:DataFrame
using ProgressMeter
import CSV:CSV
using Random
# Note: Set JULIA_NUM_THREADS to the amount of threads to use.

# Number of runs, per approach, per instance
# n_exp = 1
# (Maximum) amount of time for each run, per instance in seconds.
t_max = 1.0 # 1 s, can also be set to Inf to be unbounded.
# (Maximum) amount of evaluations

# Set to the right value!
e_max_lookup = Dict(
    # Format `(number of orders, number of machines) => number of evaluations,`
    (100, 5)  => 235879800,
    (100, 10) => 266211000,
    (100, 20) => 283040000,
    (200, 10) => 272515500,
    (200, 20) => 287728850,
    (500, 20) => 260316750,
    (20, 5)   => 182224100,
    (20, 10)  => 224784800,
    (20, 20)  => 256896400,
    (50, 5)   => 220712150,
    (50, 10)  => 256208100,
    (50, 20)  => 275954150,
)
# Sidenote: An approach can converge and not use up the evaluations.

moments = [t_max]

if length(ARGS) >= 3
    exp_idx = parse(Int64, ARGS[1])
    instance_idx = parse(Int64, ARGS[2])
    approach_idx = parse(Int64, ARGS[3])
else
    error("julia -O3 run_single_taillard_per_cat_eval_limit.jl [exp_idx] [instance_idx] [approach_idx]")
end

if length(ARGS) >= 4
    exp_count = parse(Int64, ARGS[4])
else
    exp_count = 1
end

path_results_time = "./results/results_pfs_taillard_instance$(instance_idx)_approach$(approach_idx)_exp$(exp_idx)_time.csv"
path_results_evals = "./results/results_pfs_taillard_instance$(instance_idx)_approach$(approach_idx)_exp$(exp_idx)_evals.csv"
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
    ("Random Key SimpleGA", (f, n, t, e) -> optimize_rksimplega(f, n, t, e)),
    
    # Integer Permutation SimpleGA with various permutation crossover operators.
    ("Integer Permutation SimpleGA - PMX", 
        (f, n, t, e) -> optimize_ipsimplega(PMX(n), f, n, t, e)),
    ("Integer Permutation SimpleGA - OX", 
        (f, n, t, e) -> optimize_ipsimplega(OX(n), f, n, t, e)),
    ("Integer Permutation SimpleGA - CX", 
        (f, n, t, e) -> optimize_ipsimplega(CX(n), f, n, t, e)),
    ("Integer Permutation SimpleGA - ER", 
        (f, n, t, e) -> optimize_ipsimplega(ER(n), f, n, t, e)),
]

println("Reading and parsing instances")
# Load and parse all instances.
instances = collect(Iterators.flatten(
    begin 
        instance_file = open(instance_path, "r"); 
        instance_str = read(instance_file, String);
        close(instance_file);
        ( (string(basename(instance_path), "#", i), instance) 
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
    println("Warming up approach")
    instance_warmup = instances[1][2]
    bb_warmup = bb_wrap_pfs(instance_warmup)
    let (approach_name, optimize_approach) = approaches[approach_idx]
        # ~2.0s should be sufficient to get Julia to compile and run everything.
        # Alternatively, 5000 evaluations indicate that most things have been ran as well.
        optimize_approach(bb_warmup, instance_warmup.n, 2.0, 5000)
    end
    println("Starting experiment.")
  
    sort!(moments)
    
    let (instance_name, instance) = instances[instance_idx], 
        (approach_name, optimize_approach) = approaches[approach_idx]
    for exp_i in exp_idx:(exp_idx+exp_count-1)

        e_max = e_max_lookup[instance.n, instance.m]
        moments_eval = [e_max]  
        sort!(moments_eval)

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
        for (time, fitness) in zip(trace.moments, trace.results)
            push!(results_time, (instance_name, approach_name, exp_i, time, -fitness))
        end
        for (evals, fitness) in zip(trace.moments_n_eval, trace.results_eval)
            push!(results_eval, (instance_name, approach_name, exp_i, evals, -fitness))
        end
    end
    end
    # Store results to disk.
    results_time |> CSV.write(path_results_time)
    results_eval |> CSV.write(path_results_evals)
end