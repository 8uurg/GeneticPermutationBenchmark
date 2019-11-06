println("Setting up...")
import Glob:glob
import DataFrames:DataFrame
using ProgressMeter
import CSV:CSV

# Number of runs, per approach, per instance
n_exp = 10
# (Maximum) amount of time for each run, per instance.
t_max = 10.0

path_results_time = "./results/results_time_$(ARGS[1])_$(ARGS[2]).csv"
path_instances = "./instances/qaplib/instances"

# Make sure ./results exists
if !isdir("./results")
    mkdir("./results")
end

# Find instance files.
instances = glob(string(path_instances, "/*.dat"))

# Check if instances are present.
err_no_instances = """No instances found.
Download instances using `fetch.py` in ./instances"""
@assert length(instances) != 0 err_no_instances

println("Loading problem & approaches...")
# Load problem utilities
include("./qap.jl")
# Load performance tracking utilities
include("../../utilities/trace.jl")
# Load approaches
include("../../approaches/pGOMEA/pGOMEA.jl")
include("../../approaches/qGOMEA/qGOMEA.jl")
include("../../approaches/SimpleGA/IPSimpleGA.jl")
include("../../approaches/SimpleGA/RKSimpleGA.jl")

# And set them up
approaches = [
    ("Permutation GOMEA", (f, n, t) -> optimize_pgomea(f, n, t)),
    ("qGOMEA", (f, n, t) -> optimize_qgomea(f, n, t)),
    ("Random Key SimpleGA", (f, n, t) -> optimize_rksimplega(f, n, t)),
    ("Integer Permutation (PMX) SimpleGA", 
        (f, n, t) -> optimize_ipsimplega(PMX(n), f, n, t)),
    ("Integer Permutation (OX) SimpleGA", 
        (f, n, t) -> optimize_ipsimplega(OX(n), f, n, t)),
    ("Integer Permutation (CX) SimpleGA", 
        (f, n, t) -> optimize_ipsimplega(CX(n), f, n, t)),
    ("Integer Permutation (ER) SimpleGA", 
        (f, n, t) -> optimize_ipsimplega(ER(n), f, n, t)),
]

if length(ARGS) >= 2
    approaches = approaches[parse(Int64, ARGS[1]):parse(Int64, ARGS[2])]
end

println("Reading and parsing instances")
# Load and parse all instances.
instances = [begin 
    instance_file = open(instance_path, "r"); 
    instance_str = read(instance_file, String);
    close(instance_file);
    (basename(instance_path), parse_qap_qaplib(Int64, instance_str))
end for instance_path in instances]

# Initialize results storage
results_time = DataFrame(
    [   String,    String,  Int64, Float64,    Float64], 
    [:instance, :approach, :exp_i,   :time, :objective])

moments = [0.1, 0.5, 1.0, 5.0, 10.0]
println("Starting experiment")
progress = Progress(length(instances)*length(approaches)*n_exp)
@showprogress for (instance_name, instance) in instances
    bbf = bb_wrap_qap(instance)
    for (approach_name, optimize_approach) in approaches, exp_i in 1:n_exp
        # Perform GC for good measure.
        GC.gc()
        # Setup performance/time trace in black box.
        trace = ProgressTrace(moments, 0.0)
        bbf_traced = trace_bb(bbf, trace)
        # Run experiment
        res = optimize_approach(bbf_traced, instance.n, t_max)
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