println("Setting up...")
import Glob:glob
import DataFrames:DataFrame
using ProgressMeter
import CSV:CSV

# Number of runs, per approach, per instance
n_exp = 10
# (Maximum) amount of time for each run, per instance.
t_max = Inf64 # 10.0
# (Maximum) amount of evaluations
e_max = 10000000

moments = [t_max]
moments_eval = [e_max]

path_results_time = "./results/results_time_$(ARGS[1])_$(ARGS[2]).csv"
path_instances = "./instances/taillard/instances"

# Make sure ./results exists
if !isdir("./results")
    mkdir("./results")
end

# Find instance files.
instances = glob(string(path_instances, "/*.dat"))

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

if length(ARGS) >= 2
    approaches = approaches[parse(Int64, ARGS[1]):parse(Int64, ARGS[2])]
end

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

begin
    println("Warming up approaches")
    bb_warmup = bb_wrap_pfs(instances[1])
    for (approach_name, optimize_approach) in approaches
        # ~2.0s should be sufficient to get Julia to compile and run everything.
        # Alternatively, 5000 evaluations indicate that most things have been ran as well.
        optimize_approach(bb_warmup, instances[1].n, 2.0, 5000)
    end
    println("Starting experiment")
    progress = Progress(length(instances)*length(approaches)*n_exp)
    @showprogress for (instance_name, instance) in instances
        bbf = bb_wrap_pfs(instance)
        # Test evaluation.
        bbf(shuffle!(collect(1:instance.n)))
        for (approach_name, optimize_approach) in approaches, exp_i in 1:n_exp
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
            next!(progress)
        end
        # Store results
        results_time |> CSV.write(path_results_time)
    end
end