using Random
import Statistics:quantile
include("../../utilities/FastClustering.jl")
include("../populationsizingscheme.jl")

## Inverse permutation with preset location
function wrap_rkeys_to_permutation(f :: Function)
    prm = collect(1:100)
    function ev(assignment :: Vector{Float64})
        resize!(prm, length(assignment))
        sortperm!(prm, assignment)
        return f(prm)
    end
    return ev
end

## Circular shifting of a subset through reversals.
function circshift_inplace!(A, shifts, start=firstindex(A), stop=lastindex(A))
     shifts = mod(shifts, length(A))
     reverse!(A, start, stop)
     @inbounds reverse!(A, start, start+shifts-1)
     @inbounds reverse!(A, start+shifts, stop)
     return A
end

circshift1_inplace!(A) = circshift_inplace!(A, 1)


## Recombination of position
function mix!(dst :: Vector{Float64}, src :: Vector{Float64}, mask :: Vector{Int64})
    for i in mask
        dst[i] = src[i]
    end
    return dst
end

function reencode!(sol :: Vector{Float64}, perm :: Vector{Int64}, rk :: Vector{Float64}, rng :: MersenneTwister)
    # Generate new random keys
    @inbounds for i in 1:length(rk)
        rk[i] = rand(rng, Float64)
    end
    sort!(rk)
    # Find the permutation that sorts the original keys.
    sortperm!(perm, sol)
    # Inversely permute the new keys in rk onto sol.
    # So that the original `sol` and new `sol` both encode the same underlying permutation.
    @inbounds for (i, j) in enumerate(perm)
        # The key at index j has rank i. `rk` is sorted, thus index = rank.
        sol[j] = rk[i] 
    end
    return sol
end

function random_rescale!(sol :: Vector{Float64}, mask :: Vector{Int64}, rng :: MersenneTwister)
    # Find the minimum and maximum to rescale
    sm_min, sm_max = extrema(sol[i] for i in mask)
    sm_range = sm_max - sm_min
    # How many intervals?
    n_intervals = length(sol)
    # Pick an interval [random_left_bracket:random_left_bracket+1/n_intervals]
    # to scale elements to.
    random_left_bracket = rand(rng, 0:(n_intervals-1)) / n_intervals
    # NOTE! In order to avoid equal keys, and hence potential bias towards the solution [1, ..., n]
    # Use a random float rather than an integer.
    # random_left_bracket = rand(rng, Float64) * (n_intervals - 1) / n_intervals
    
    # Remap keys from [sm_min, sm_max] to [random_left_bracket:random_left_bracket+1/n_intervals]
    if sm_range == 0
        for i in mask
            sol[i] = random_left_bracket
        end
    else
        for i in mask
            sol[i] = ((sol[i] - sm_min) / sm_range) * (1.0 / n_intervals) + random_left_bracket
        end
    end

    return sol
end

##
mutable struct PGomeaSolution
    perm :: Vector{Float64}
    fitness :: Float64
end
Base.isless(a :: PGomeaSolution, b :: PGomeaSolution) = isless(a.fitness, b.fitness)
Base.isequal(a :: PGomeaSolution, b :: PGomeaSolution) = a.perm == b.perm
##
function population_has_converged(pop :: Vector{X}) :: Bool where {X}
    if length(pop) <= 1
        return true
    end
    frst = first(pop)
    return all(isequal(frst, s) for s in Iterators.drop(pop, 1))
end

##
struct PGomeaMixer
    # Important
    f :: Function
    n :: Int64
    population :: Vector{PGomeaSolution}
    D :: Matrix{Float64}
    fos :: Vector{Vector{Int64}}
    # Internal config.
    forced_improvement :: Symbol
    fos_type :: Union{Symbol, Vector{Vector{Int64}}}
    rescale_probability :: Float64
    crf :: ClusteringReductionFormula
    # Stats & info
    best :: Ref{PGomeaSolution}

    generations :: Ref{Int64}
    generations_no_improvement :: Ref{Int64}
    converged :: Ref{Bool}
    # Memory allocations.
    mixing_backup :: Vector{Float64}
    
    reencode_keys :: Vector{Float64}
    reencode_perm :: Vector{Int64}

    # RNG
    rng :: MersenneTwister

    function PGomeaMixer(
        f :: Function,
        n :: Int64,
        population :: Vector{PGomeaSolution},
        D :: Matrix{Float64},
        fos :: Vector{Vector{Int64}},
        forced_improvement :: Symbol,
        fos_type :: Union{Symbol, Vector{Vector{Int64}}},
        rescale_probability :: Float64,
        crf :: ClusteringReductionFormula,
        rng :: MersenneTwister)
        #
        new(f, n, population, D, fos,
            forced_improvement, fos_type, rescale_probability, crf,
            Ref(maximum(population)),
            Ref(0), Ref(0), Ref(false),
            collect(1:n), 
            collect(LinRange(0.0, 1.0, n)), collect(1:n),
            rng)
    end

    function PGomeaMixer(
        f :: Function,
        n :: Int64,
        population :: Vector{PGomeaSolution},
        D :: Matrix{Float64},
        fos :: Vector{Vector{Int64}},
        forced_improvement :: Symbol,
        fos_type :: Union{Symbol, Vector{Vector{Int64}}},
        rescale_probability :: Float64,
        crf :: ClusteringReductionFormula,
        best :: Ref{PGomeaSolution},
        rng :: MersenneTwister)
        #
        best[] = max(best[], maximum(population))
        new(f, n, population, D, fos,
            forced_improvement, fos_type, rescale_probability, crf,
            best,
            Ref(0), Ref(0), Ref(false),
            collect(1:n), 
            collect(LinRange(0.0, 1.0, n)), collect(1:n),
            rng)
    end
end

function evaluate(f :: Function, perm :: Vector{Int64})
    fitness = f(perm)
    return PGomeaSolution(perm, fitness)
end

function calcD!(pm :: PGomeaMixer)
    fill!(pm.D, zero(Float64))
    @inbounds for individual in pm.population
        for i in 1:pm.n
            for j in i+1:pm.n
                pm.D[i, j] += abs(individual.perm[i] - individual.perm[j])
                # pm.D[i, j] += (individual.perm[i] - individual.perm[j])^2
                pm.D[j, i] = pm.D[i, j]
            end
        end
    end
    #pm.D .= maximum(pm.D) .- pm.D
    pm.D
end

function entropy(x) 
    if x > 0.0 && x < 1.0 
        return -x * log2(x) + - (1.0-x) * log2(1.0-x)
    else 
        return 0.0
    end
end

double_entropy(x) = entropy(entropy(x))

weighted_once_inv_double_entropy(x, w=0.5) = w*entropy(x)+(1.0-w)*(1.0-entropy(entropy(x)))
inv_weighted_once_inv_double_entropy(x, w=0.5) = 1.0 - weighted_once_inv_double_entropy(x, w)

function calcD_original_pw!(pm :: PGomeaMixer)
    fill!(pm.D, zero(Float64))
    @inbounds for i in 1:pm.n
        for j in i+1:pm.n
            c_ab = 0
            c_dist = 0
            for individual in pm.population
                c_dist += abs(individual.perm[i] - individual.perm[j])
                c_ab += ifelse(individual.perm[i] > individual.perm[j], 1, 0)
            end
            δ₁ = c_dist / (pm.n^2)
            δ₂ = entropy(c_ab / pm.n)
            pm.D[i, j] = δ₁ * δ₂
            pm.D[j, i] = pm.D[i, j]
        end
    end
    #pm.D .= maximum(pm.D) .- pm.D
    pm.D
end

function calcD_original!(pm :: PGomeaMixer)
    fill!(pm.D, zero(Float64))
    @inbounds for i in 1:pm.n
        for j in i+1:pm.n
            c_ab = 0
            c_dist = 0
            for individual in pm.population
                c_dist += (individual.perm[i] - individual.perm[j])^2
                c_ab += ifelse(individual.perm[i] > individual.perm[j], 1, 0)
            end
            δ₁ = 1 - (sqrt(c_dist) / pm.n)
            δ₂ = 1 - entropy(c_ab / pm.n)
            # Invert direction.
            pm.D[i, j] = -1 * δ₁ * δ₂
            pm.D[j, i] = pm.D[i, j]
        end
    end
    #pm.D .= maximum(pm.D) .- pm.D
    pm.D
end

function calcD_original_regularized!(pm :: PGomeaMixer)
    fill!(pm.D, zero(Float64))
    @inbounds for i in 1:pm.n
        for j in i+1:pm.n
            c_ab = 0
            c_dist = 0
            for individual in pm.population
                c_dist += abs(individual.perm[i] - individual.perm[j])
                c_ab += ifelse(individual.perm[i] > individual.perm[j], 1, 0)
            end
            δ₁ = c_dist / (pm.n^2)
            δ₂ = inv_weighted_once_inv_double_entropy(c_ab / pm.n, 0.3) + 0.3
            pm.D[i, j] = δ₁ * δ₂
            pm.D[j, i] = pm.D[i, j]
        end
    end
    #pm.D .= maximum(pm.D) .- pm.D
    pm.D
end

function calcD_random!(pm :: PGomeaMixer)
    for i in 1:pm.n
        for j in i+1:pm.n
            pm.D[i, j] = rand(pm.rng)
            pm.D[j, i] = pm.D[i, j]
        end
    end
    pm.D
end

function edamixing(sol :: PGomeaSolution, pm :: PGomeaMixer; shuffle_fos=true, donor_fixed :: Union{Nothing, Ref{PGomeaSolution}} = nothing)
    # Shuffle if required.
    if shuffle_fos
        shuffle!(pm.rng, pm.fos)
    end
    # Alias for convinience.
    dst = sol.perm
    copyto!(pm.mixing_backup, dst)
    # Keep track of some informative stats.
    solution_changed = false
    current_improved = false
    best_improved = false
    # Actual mixing.
    for s in pm.fos
        if donor_fixed === nothing
            donor_sol = rand(pm.rng, pm.population)
        else
            donor_sol = donor_fixed[]
        end
        donor = donor_sol.perm
        
        # If subsets are equal in value, no mixing is going to make
        # a difference
        if all(pm.mixing_backup[i] == donor[i] for i in s)
            continue
        end

        mix!(pm.mixing_backup, donor, s)

        if rand(pm.rng) < pm.rescale_probability
            random_rescale!(pm.mixing_backup, s, pm.rng)
        end
        
        if pm.mixing_backup == dst
            # Solution stayed the same... No reevaluation needed.
            continue
        end
        fitness = pm.f(pm.mixing_backup)
        if fitness < sol.fitness #|| (pm.mixing_backup == donor)
            copyto!(pm.mixing_backup, dst)
        elseif fitness == sol.fitness
            # GOMEA generally copies over even without improvements.
            copyto!(dst, pm.mixing_backup)
            if dst != pm.mixing_backup
                solution_changed = true
            end
        elseif fitness > sol.fitness
            copyto!(dst, pm.mixing_backup)
            sol.fitness = fitness
            current_improved = true
            solution_changed = true
            if fitness > pm.best[].fitness
                pm.best[] = sol
                best_improved = true
            end
        end
    end
    return best_improved, solution_changed
end

function step!(pm :: PGomeaMixer)
    # Re-encode population
    for p in pm.population
        reencode!(p.perm, pm.reencode_perm, pm.reencode_keys, pm.rng)
    end

    # Calculate D -- the Dependency Matrix -- using the population.
    is_linkage_tree = false
    is_swapping = false
    if pm.fos_type == :distance
        calcD!(pm)
        is_linkage_tree = true
    elseif pm.fos_type == :distance_swap
        calcD!(pm)
        is_linkage_tree = true
        is_swapping = true
    elseif pm.fos_type == :original
        calcD_original!(pm)
        is_linkage_tree = true
    elseif pm.fos_type == :original_swap
        calcD_original!(pm)
        is_linkage_tree = true
        is_swapping = true
    elseif pm.fos_type == :random
        calcD_random!(pm)
        is_linkage_tree = true
    elseif pm.fos_type == :univariate
        # Univariate Factorization
        empty!(pm.fos)
        append!(pm.fos, [i] for i in 1:pm.n)
    elseif typeof(pm.fos_type) <: Vector{Vector{Int64}}
        # Predetermined Factorization
        empty!(pm.fos)
        append!(pm.fos, pm.fos_type)
    else
        error("Unknown FOS type $(pm.fos_type)")
    end
    # calcD_original_regularized!(pm)
    
    # Compute FOS using D
    if is_linkage_tree
        empty!(pm.fos)
        parent_idx = zeros(Int64, 2*pm.n-1)
        #fos_indexset = LCP(pm.D, pm.crf, pm.rng; parent_idx=parent_idx)
        fos_indexset = LCP(pm.D, pm.crf, pm.rng; parent_idx=parent_idx, randomized=true)
        append!(pm.fos, collect(a) for (i,a) in enumerate(fos_indexset))
        if is_swapping
            rewrite_by_swap_fos(pm.fos, findchildrenpairs(parent_idx), pm.D)
        end
    end

    improved_any = false

    # Filter the FOS.
    # sort!(pm.fos, by=f->length(f))
    # filter!(f -> (length(f) > 1 && length(f) < round(pm.n/12*10)), pm.fos)
    # filter!(f -> (length(f) < round(pm.n/12*10)), pm.fos)
    # filter!(f -> (length(f) > 1), pm.fos)
    # filter!(f -> (length(f) > 1 && length(f) < pm.n), pm.fos)
    filter!(f -> (length(f) < pm.n), pm.fos)
    
    # Determine Forced Improvement threshold.
    fi_threshold = typemax(Int64) # :)
    allow_fi_upon_unchanged_solution = false
    if pm.forced_improvement == :original || pm.forced_improvement == :original_sc
        fi_threshold = floor(Int64, 1+log10(length(pm.population)))
    elseif pm.forced_improvement == :extended || pm.forced_improvement == :extended_sc
        fi_threshold = 10 + 10*floor(Int64, log10(length(pm.population)))
    elseif pm.forced_improvement == :basically_none || pm.forced_improvement == :basically_none_sc
        fi_threshold = 10 + pm.n*floor(Int64, log10(length(pm.population)))
    elseif pm.forced_improvement == :none || pm.forced_improvement == :none_sc
    else
        error("Invalid FI strategy.")
    end
    if pm.forced_improvement == :original_sc || pm.forced_improvement == :extended_sc ||
        pm.forced_improvement == :basically_none_sc || pm.forced_improvement == :none_sc
        allow_fi_upon_unchanged_solution = true
    end

    # Perform Optimal Mixing on population
    for individual in pm.population
        # Normal Mixing.
        improved_any_this_mix, solution_changed = edamixing(individual, pm)
        improved_any |= improved_any_this_mix
        # Forced improvement.
        if pm.forced_improvement != :none && ((allow_fi_upon_unchanged_solution && !solution_changed) ||
            pm.generations_no_improvement[] > fi_threshold)
            improved_any_this_mix, _ = edamixing(individual, pm; donor_fixed = pm.best)
            improved_any |= improved_any_this_mix
        end
    end

    # Update statistics
    pm.generations[] += 1

    if improved_any
        pm.generations_no_improvement[] = 0
    else
        pm.generations_no_improvement[] += 1
    end

    if population_has_converged(pm.population) || pm.generations_no_improvement[] >  1 + 10*floor(Int64, log10(length(pm.population)))*20
        pm.converged[] = true
        # println("A population of size $(length(pm.population)) has converged.")
    end

    return improved_any
end

function generate_new_pgomeasolution_random(f :: Function, n :: Int64, rng :: MersenneTwister)
    perm = rand(rng, Float64, n)
    PGomeaSolution(perm, f(perm))
end

function generate_new_pgomeasolution_bounds(bounds :: Vector{Tuple{Float64, Float64}})
    function generate_new_solution(f :: Function, n :: Int64, rng :: MersenneTwister)
        @assert length(bounds) == n
        perm = [rand(rng) * (a[2]-a[1]) + a[1] for a in bounds]
        PGomeaSolution(perm, f(perm))
    end
    return generate_new_solution
end

function generate_new_pgomeasolution_bounds(bounds :: Tuple{Vector{Float64}, Vector{Float64}})
    function generate_new_solution(f :: Function, n :: Int64, rng :: MersenneTwister)
        @assert length(bounds) == n
        perm = invperm(sortperm([rand(rng) * (ub-lb) + lb for (lb, ub) in zip(bounds[1], bounds[2])]))
        PGomeaSolution(perm, f(perm))
    end
    return generate_new_solution
end

function create_mixer(f :: Function, n :: Int64, population_size :: Int64, 
    forced_improvement :: Symbol,
    fos_type :: Union{Symbol, Vector{Vector{Int64}}},
    rescale_probability :: Float64,
    crf :: ClusteringReductionFormula,
    rng :: MersenneTwister,
    best :: Union{Nothing, Ref{PGomeaSolution}} = nothing;
    initial_solution_generator :: Function) :: PGomeaMixer
    # Generate initial population.
    population = [initial_solution_generator(f, n, rng) for _ in 1:population_size]
    # Create matrices and FoS vector.
    D = zeros(Float64, (n, n))
    fos = Vector{Vector{Int64}}()
    sizehint!(fos, 2*n)
    if best === nothing
        return PGomeaMixer(f, n, population, D, fos, forced_improvement, fos_type, rescale_probability, crf, rng)
    else
        return PGomeaMixer(f, n, population, D, fos, forced_improvement, fos_type, rescale_probability, crf, best, rng)
    end
end

function optimize_pgomea(rf :: Function, n :: Int64, t=10.0, e=typemax(Int64);
    initial_solution_generator :: Function = generate_new_pgomeasolution_random,
    population_size_base=4, population_sizing_step_factor :: Int64 = 4,
    rescale_probability :: Float64 = 0.1,
    crf=UPGMA(),
    forced_improvement :: Symbol = :extended,
    fos_type :: Union{Symbol, Vector{Vector{Int64}}} = :original,
    target_fitness :: Union{Nothing, Float64} = nothing)
    #
    time_start = time()
    n_evals = 0

    rng = MersenneTwister()
    
    fx = wrap_rkeys_to_permutation(rf)
    function f(sol :: Vector{Float64})
        n_evals += 1
        return fx(sol)
    end

    next_population_size = population_size_base*2

    initial_mixer = create_mixer(f, n, population_size_base, forced_improvement, fos_type, rescale_probability, crf, rng, initial_solution_generator=initial_solution_generator)
    mixers = PGomeaMixer[initial_mixer]
    steps = 0
    last_steps = 0
    best = initial_mixer.best
    upto_gen = SubpopulationStepcountGenerator(population_sizing_step_factor)

    while (time() - time_start < t) && (n_evals <= e) && best[].fitness != target_fitness
        upto_mixer, _ = iterate(upto_gen)
        upto_mixer = min(upto_mixer, length(mixers))

        for i_mixer in 1:upto_mixer
            if i_mixer == length(mixers)
                last_steps += 1
            end

            step!(mixers[i_mixer])
        end

        # Add new metaheuristic upon having stepped the last one 4 times.
        # Or if all have converged. Oh well.
        if last_steps == population_sizing_step_factor || length(mixers) == 0
            last_steps = 0
            push!(mixers, create_mixer(f, n, next_population_size, forced_improvement, fos_type, rescale_probability, crf, rng, best, initial_solution_generator=initial_solution_generator))
            step!(mixers[end])
            # println("Created new population of size $(next_population_size)")
            next_population_size *= 2
        end
        filter!(f -> !f.converged[], mixers)

        steps += 1
    end

    return (best[].fitness, sortperm(best[].perm))
end
