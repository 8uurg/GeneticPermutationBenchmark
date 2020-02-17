println("Setting up...")
using Distributed
n_proc = 8
addprocs(n_proc)
@everywhere import Glob:glob
@everywhere import DataFrames:DataFrame
@everywhere using ProgressMeter
@everywhere import CSV:CSV
@everywhere using Random

# Note: Set JULIA_NUM_THREADS to the amount of threads to use.

# Number of runs, per approach, per instance
@everywhere n_exp = 25
@everywhere success_threshold = 5
# (Maximum) amount of time for each run, per instance in seconds.
@everywhere t_max = 100.0
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
path_results_time = "./results/results_progressive_exphit_biv_invonly_ord_$(exp_idx_offset)_time.csv"
path_results_evals = "./results/results_progressive_exphit_biv_invonly_ord_$(exp_idx_offset)_evals.csv"
path_results_time_hit = "./results/results_progressive_exphit_biv_invonly_ord_$(exp_idx_offset)_time_hit.csv"
path_results_evals_hit = "./results/results_progressive_exphit_biv_invonly_ord_$(exp_idx_offset)_evals_hit.csv"

# Make sure ./results exists
if !isdir("./results")
    mkdir("./results")
end

println("Loading problem & approaches...")
# Load problem utilities
@everywhere include("./biv.jl")
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
    ("qGOMEA - LT/Distance - No FI - OX", 
        (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, forced_improvement = :none, target_fitness=target_fitness)),
    
    ("qGOMEA - LT/PermutationGOMEA Original", 
        (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, fos_type=:original, target_fitness=target_fitness)),
    # ("qGOMEA - LT/Distance - 10x FI - PMX", 
    #     (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, permutation_repair=:pmx, target_fitness=target_fitness)),

    # ("qGOMEA - RT - 1x FI - OX", 
    #     (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, forced_improvement = :original, fos_type=:random, target_fitness=target_fitness)),
    ("qGOMEA - RT - 10x FI - OX", 
        (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, forced_improvement = :extended, fos_type=:random, target_fitness=target_fitness)),
    ("qGOMEA - RT - No FI - OX", 
        (f, n, t, e; target_fitness) -> optimize_qgomea(f, n, t, e, forced_improvement = :none, fos_type=:random, target_fitness=target_fitness)),
    
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
@everywhere func = [
    ("Inversion", sorted_inversion, n -> convert(Float64, div( n*( n-1), 2)) ), 
    # ("Sequential Inversion", sorted_sequential_inversion, n -> convert(Float64, n - 1)),
    # ("Sequential Pairs", sorted_sequential_pairs, n -> convert(Float64, n - 1))
]
@everywhere ns = [10, 15, 20, 25, 50, 100, 200, 400, 800]

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
    @everywhere instance_warmup = BIVInstance(sorted_inversion, 10)
    @everywhere bb_warmup = bb_wrap_biv(instance_warmup)
    @everywhere for (approach_name, optimize_approach) in approaches
        # ~2.0s should be sufficient to get Julia to compile and run everything.
        # Alternatively, 5000 evaluations indicate that most things have been ran as well.
        optimize_approach(bb_warmup, instance_warmup.n, 2.0, 1000; target_fitness=nothing)
    end
    @everywhere instance_warmup = BIVInstance(sorted_sequential_pairs, 10)
    @everywhere bb_warmup = bb_wrap_biv(instance_warmup)
    @everywhere for (approach_name, optimize_approach) in approaches
        # ~2.0s should be sufficient to get Julia to compile and run everything.
        # Alternatively, 5000 evaluations indicate that most things have been ran as well.
        optimize_approach(bb_warmup, instance_warmup.n, 2.0, 1000; target_fitness=nothing)
    end
    @everywhere instance_warmup = BIVInstance(sorted_sequential_inversion, 10)
    @everywhere bb_warmup = bb_wrap_biv(instance_warmup)
    @everywhere for (approach_name, optimize_approach) in approaches
        # ~2.0s should be sufficient to get Julia to compile and run everything.
        # Alternatively, 5000 evaluations indicate that most things have been ran as well.
        optimize_approach(bb_warmup, instance_warmup.n, 2.0, 1000; target_fitness=nothing)
    end
    
    println("Getting ready.")
    @everywhere expected_exps = length(approaches) * length(func) * length(ns) * n_exp
    lc_progress_channel = RemoteChannel(()->Channel{Int}(expected_exps))
    @everywhere progress_channel = $lc_progress_channel

    @everywhere function run_singular_exp(instance :: BIVInstance, instance_opt :: Float64, approach_idx :: Int64)
        optimize_approach = approaches[approach_idx][2]
        bbf = bb_wrap_biv(instance)
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
        put!(progress_channel, 1)
        return trace, res
    end

    @everywhere function run_experiment_tier(func_idx :: Int64, n_idx :: Int64, approach_idx :: Int64)
        benchmark_func_name, benchmark_func, benchmark_func_n_to_opt = func[func_idx]
        n = ns[n_idx]
        instance_name = "$benchmark_func_name n=$n"
        # println("Running $instance_name for approach $approach_idx")
        instance = BIVInstance(benchmark_func, n)
        instance_opt = benchmark_func_n_to_opt(n)
        exps = fetch.([@spawnat :any run_singular_exp(instance, instance_opt, approach_idx) for _ in 1:n_exp])
        successes = 0
        # println(exps)
        for (trace, result) in exps
            # Do some tier based processing
            if trace.hitting_time < Inf
                successes += 1
            end
        end
        # println("$instance_name for $approach_idx with $successes successes.")
        if successes < success_threshold
            # Send notification
            put!(progress_channel, n_exp * (length(ns) - n_idx + 1) - n_exp)
            # Too many failures, return without further work.
            return [(instance_name, instance_opt, exps)]
        end

        exps_more_items = Future[]
        if n_idx < length(ns)
            exps_more_items = @spawnat :any run_experiment_tier(func_idx, n_idx + 1, approach_idx)
        end

        return reduce(vcat, [[(instance_name, instance_opt, exps)], fetch(exps_more_items)]) 
    end
        
    
    println("Starting experiment, running a maximum of $(expected_exps) experiments on $(nworkers()) worker(s).")

    progress = Progress(expected_exps)

    local exps_approaches :: Any
    
    @sync begin
        @async begin
            exps_approaches_futures = [(i_a, @spawnat :any run_experiment_tier(i_f, 1, i_a)) for i_a in 1:length(approaches) for i_f in 1:length(func)]
            exps_approaches = map(f -> (f[1], fetch(f[2])), exps_approaches_futures)
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
    for (approach_idx, exps_per_instance) in exps_approaches
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
