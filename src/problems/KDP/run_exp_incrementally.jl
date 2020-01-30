println("Setting up...")
using Distributed
n_proc = 1
addprocs(n_proc)
@everywhere import Glob:glob
@everywhere import DataFrames:DataFrame
@everywhere using ProgressMeter
@everywhere import CSV:CSV
@everywhere using Random

# Note: Set JULIA_NUM_THREADS to the amount of threads to use.

# Number of runs, per approach, per instance
@everywhere n_exp = 1
@everywhere success_threshold = 1 
# (Maximum) amount of time for each run, per instance in seconds.
@everywhere t_max = 10.0
# (Maximum) amount of evaluations
@everywhere e_max = typemax(Int64) # 10000000
# Sidenote: An approach can converge and not use up the evaluations.

@everywhere moments = [t_max]
@everywhere moments_eval = [e_max]

if length(ARGS) > 0
    exp_idx_offset = parse(Int64, ARGS[1])
else
    exp_idx_offset = 0
end
path_results_time = "./results/results_progressive_exphit_kdp_loose_shuf_$(exp_idx_offset)_time.csv"
path_results_evals = "./results/results_progressive_exphit_kdp_loose_shuf_$(exp_idx_offset)_evals.csv"
path_results_time_hit = "./results/results_progressive_exphit_kdp_loose_shuf_$(exp_idx_offset)_time_hit.csv"
path_results_evals_hit = "./results/results_progressive_exphit_kdp_loose_shuf_$(exp_idx_offset)_evals_hit.csv"


# Make sure ./results exists
if !isdir("./results")
    mkdir("./results")
end

println("Loading problem & approaches...")
# Load problem utilities
@everywhere include("./kdp.jl")
# Load performance tracking utilities
@everywhere include("../../utilities/trace.jl")
# Include permutation remapper
@everywhere include("../permutationtools.jl")
# Load approaches
@everywhere include("../../approaches/pGOMEA/pGOMEA.jl")
@everywhere include("../../approaches/qGOMEA/qGOMEA.jl")
@everywhere include("../../approaches/SimpleGA/IPSimpleGA.jl")
@everywhere include("../../approaches/SimpleGA/RKSimpleGA.jl")

# And set them up
@everywhere approaches = [
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
    # ("qGOMEA - LT/Distance - No FI - OX", 
    #     (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, forced_improvement = :none, target_fitness=target_fitness)),
    
    ("qGOMEA - LT/PermutationGOMEA Original", 
        (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, fos_type=:original, target_fitness=target_fitness)),
    # ("qGOMEA - LT/Distance - 10x FI - PMX", 
    #     (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, permutation_repair=:pmx, target_fitness=target_fitness)),

    # ("qGOMEA - RT - 1x FI - OX", 
    #     (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, forced_improvement = :original, fos_type=:random, target_fitness=target_fitness)),
    ("qGOMEA - RT - 10x FI - OX", 
        (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, forced_improvement = :extended, fos_type=:random, target_fitness=target_fitness)),
    # ("qGOMEA - RT - No FI - OX", 
    #     (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, forced_improvement = :none, fos_type=:random, target_fitness=target_fitness)),
    
    # Random Key SimpleGA
    ("Random Key SimpleGA", (f, n, t, e; target_fitness) -> optimize_rksimplega(f, n, t, e, target_fitness=target_fitness)),
    
    # Integer Permutation SimpleGA with various permutation crossover operators.
    ("Integer Permutation SimpleGA - PMX", 
        (f, n, t, e; target_fitness) -> optimize_ipsimplega(PMX(n), f, n, t, e, target_fitness=target_fitness)),
    ("Integer Permutation SimpleGA - OX", 
        (f, n, t, e; target_fitness) -> optimize_ipsimplega(OX(n), f, n, t, e, target_fitness=target_fitness)),
    ("Integer Permutation SimpleGA - CX", 
        (f, n, t, e; target_fitness) -> optimize_ipsimplega(CX(n), f, n, t, e, target_fitness=target_fitness)),
    ("Integer Permutation SimpleGA - ER", 
        (f, n, t, e; target_fitness) -> optimize_ipsimplega(ER(n), f, n, t, e, target_fitness=target_fitness)),
]

println("Initializing instances")
# Initialize instances.
@everywhere n_blocks = [3, 4, 5, 10, 20]
@everywhere n_block_size = [4, 5, 6, 7, 8, 9]

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
    @everywhere instance_warmup = KDPInstance(3, 3)
    @everywhere bb_warmup = bb_wrap_kdp(instance_warmup)
    @everywhere for (approach_name, optimize_approach) in approaches
        # ~2.0s should be sufficient to get Julia to compile and run everything.
        # Alternatively, 5000 evaluations indicate that most things have been ran as well.
        optimize_approach(bb_warmup, instance_warmup.n, 2.0, 5000; target_fitness=nothing)
    end
    
    println("Getting ready.")
    @everywhere expected_exps = length(approaches) * length(n_blocks) * length(n_block_size) * n_exp
    lc_progress_channel = RemoteChannel(()->Channel{Int}(expected_exps))
    @everywhere progress_channel = $lc_progress_channel

    @everywhere function run_singular_exp(instance :: KDPInstance, approach_idx :: Int64)
        optimize_approach = approaches[approach_idx][2]
        instance_opt = convert(Float64, instance.b)
        bbf = random_remap(bb_wrap_kdp(instance), instance.n)
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
        put!(progress_channel, 1)
        return trace, res
    end

    @everywhere function run_experiment_tier(block_size_idx :: Int64, block_count_idx :: Int64, approach_idx :: Int64)
        block_size = n_block_size[block_size_idx]
        block_count = n_blocks[block_count_idx]
        instance_name = "k$(block_size)xb$(block_count)"
        # println("Running $instance_name for approach $approach_idx")
        instance = KDPInstance(block_size, block_count)
        instance_opt = convert(Float64, instance.b)
        exps = fetch.([@spawnat :any run_singular_exp(instance, approach_idx)])
        successes = 0
        # println(exps)
        for (trace, result) in exps
            # Do some tier based processing
            if trace.hitting_time < Inf
                successes += 1
            end
        end
        if successes < success_threshold
            # println("Not enough successes: $instance_name for $approach_idx with $successes successes.")
            # Send notification
            if block_count_idx == 1
                # We are not progressing with block size either!
                put!(progress_channel, n_exp * (length(n_blocks) - block_count_idx + 1) * (length(n_block_size) - block_size_idx + 1) - n_exp)
            else
                put!(progress_channel, n_exp * (length(n_blocks) - block_count_idx + 1) - n_exp)
            end
            # Too many failures, return without further work.
            return [(instance_name, instance_opt, exps)]
        end

        exps_more_blocks_future = Future[]
        if block_count_idx < length(n_blocks)
            # println("Requesting run with $(block_count_idx + 1) idx blocks.")
            exps_more_blocks_future = @spawnat :any run_experiment_tier(block_size_idx, block_count_idx + 1, approach_idx)
        end

        exps_larger_blocks_future = Future[]
        if block_count_idx == 1 && block_size_idx < length(n_block_size)
            # println("Requesting run with $(block_size_idx + 1) idx size blocks.")
            exps_larger_blocks_future = @spawnat :any run_experiment_tier(block_size_idx + 1, block_count_idx, approach_idx)
        end

        return reduce(vcat, [[(instance_name, instance_opt, exps)], fetch(exps_more_blocks_future), fetch(exps_larger_blocks_future)]) 
    end
        
    
    println("Starting experiment, running a maximum of $(expected_exps) experiments on $(nworkers()) worker(s).")

    progress = Progress(expected_exps)

    local exps_approaches :: Any
    
    @sync begin
        @async begin
            exps_approaches = fetch.([@spawnat :any run_experiment_tier(1, 1, i) for i in 1:length(approaches)])
            # Finish up the progressbar.
            put!(progress_channel, expected_exps)
        end
        # Progressbar Waiter
        @async begin
            local exps_done :: Int64 = 0
            while exps_done < expected_exps
                r = take!(progress_channel)
                if r == expected_exps
                    exps_done = expected_exps
                else
                    exps_done = exps_done + r
                end
                update!(progress, exps_done)
            end
        end
    end
    
    
    println("Formatting and writing results...")
    for (approach_idx, exps_per_instance) in enumerate(exps_approaches)
        approach_name = approaches[approach_idx][1]
        for (instance_name, instance_opt, exp_data) in exps_per_instance
            for (exp_i, (trace, _)) in enumerate(exp_data)
                for (time, fitness) in zip(trace.moments, trace.results)
                    push!(results_time, (instance_name, approach_name, exp_i + exp_idx_offset, time, fitness))
                end
                for (evals, fitness) in zip(trace.moments_n_eval, trace.results_eval)
                    push!(results_eval, (instance_name, approach_name, exp_i + exp_idx_offset, evals, fitness))
                end
                push!(results_time_hit, (instance_name, approach_name, exp_i + exp_idx_offset, trace.hitting_time, instance_opt))
                push!(results_eval_hit, (instance_name, approach_name, exp_i + exp_idx_offset, trace.hitting_eval, instance_opt))
            end
        end
    end

    # Store results
    results_time |> CSV.write(path_results_time)
    results_eval |> CSV.write(path_results_evals)
    results_time_hit |> CSV.write(path_results_time_hit)
    results_eval_hit |> CSV.write(path_results_evals_hit)
end
