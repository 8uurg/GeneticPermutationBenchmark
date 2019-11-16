import Random:shuffle!
import Statistics:quantile
include("../../utilities/FastClustering.jl")

## Inverse permutation with preset location
function invperm!(dst :: Vector{Int64}, src :: Vector{Int64})
    @inbounds for (i, j) in enumerate(src)
        dst[j] = i
    end
    dst
end

function invperm_ignorezero!(dst :: Vector{Int64}, src :: Vector{Int64})
    @inbounds for (i, j) in enumerate(src)
        if j != 0
            dst[j] = i
        end
    end
    dst
end

function invperm!(dst :: Vector{Int64}, src :: Vector{Int64}, subset :: Vector{Int64})
    @inbounds for i in subset
        dst[src[i]] = i
    end
    dst
end

function wrap_assignment_to_permutation(f :: Function)
    prm = collect(1:100)
    function ev(assignment :: Vector{Int64})
        resize!(prm, length(assignment))
        invperm!(prm, assignment)
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
function mix_position_ox!(dst :: Vector{Int64}, donor :: Vector{Int64}, mask :: Vector{Int64},
    storage1 :: Vector{Int64}, storage2 :: Vector{Int64})
    # ...
    fill!(storage1, -1)
    fill!(storage2, 0)
    invperm!(storage1, donor, mask)
    @inbounds dst[mask] .= 0
    invperm_ignorezero!(storage2, dst)
    q = 1
    @inbounds for j in storage2
        if j == 0
            # Skip masked items.
            continue
        end
        while storage1[q] != -1
            q += 1
        end
        storage1[q] = j
    end
    invperm!(dst, storage1)
end

function mix_position_pmx!(dst :: Vector{Int64}, donor :: Vector{Int64}, mask :: Vector{Int64},
    storage1 :: Vector{Int64}, storage2 :: Vector{Int64})
    # ...
    invperm!(storage1, dst)
    for i in mask
        odst = dst[i]
        dst[i], dst[storage1[donor[i]]] = donor[i], dst[i]
        storage1[odst], storage1[donor[i]] = storage1[donor[i]], storage1[odst]
    end
end

function create_differential_donor!(dst :: Vector{Int64}, orig_dst :: Vector{Int64}, orig_donor :: Vector{Int64}, mask :: Vector{Int64}, rng :: MersenneTwister; offset :: Union{Int64, Symbol} = 0)
    #reference = minimum((orig_donor[i],i) for i in mask)[2]
    reference = rand(rng, mask)
    if typeof(offset) === Int64
        offsett = offset
    else
        if offset == :rand_neg1to1
            offsett = rand(rng, -1:1)
        elseif offset == :rand_difference
            a = orig_dst[reference] - orig_donor[reference]
            b = 0
            a, b = minmax(a, b)
            offsett = rand(rng, a:b)*rand(rng, 0:1)
        end
    end
    for i in mask
        dst[i] = orig_donor[reference] + orig_dst[i] - orig_dst[reference] + offsett
    end
    dst_min, dst_max = extrema(dst[i] for i in mask)
    if dst_min < 1
        # Need to do a repair for negative numbers.
        for i in mask
            dst[i] += -dst_min + 1
        end
    elseif dst_max > length(dst)
        # Fix the out of bounds elements!
        dstdiff = length(dst) - dst_max
        for i in mask
            dst[i] += dstdiff
        end
    end
    return dst
end

##
mutable struct QGomeaSolution
    perm :: Vector{Int64}
    fitness :: Float64
end
Base.isless(a :: QGomeaSolution, b :: QGomeaSolution) = isless(a.fitness, b.fitness)
Base.isequal(a :: QGomeaSolution, b :: QGomeaSolution) = a.perm == b.perm
##
function population_has_converged(pop :: Vector{X}) :: Bool where {X}
    if length(pop) <= 1
        return true
    end
    frst = first(pop)
    return all(isequal(frst, s) for s in Iterators.drop(pop, 1))
end

##
struct QGomeaMixer
    # Important
    f :: Function
    n :: Int64
    population :: Vector{QGomeaSolution}
    D :: Matrix{Float64}
    fos :: Vector{Vector{Int64}}
    # Internal config.
    forced_improvement :: Symbol
    fos_type :: Symbol
    crf :: ClusteringReductionFormula
    repair :: Symbol
    # Stats & info
    best :: Ref{QGomeaSolution}

    generations :: Ref{Int64}
    generations_no_improvement :: Ref{Int64}
    converged :: Ref{Bool}
    # Memory allocations.
    mixing_backup :: Vector{Int64}
    mixing_virtual_donor :: Vector{Int64}
    mixing_perm :: Vector{Int64}
    mixing_perm_2 :: Vector{Int64}

    bs1 :: IndexSet
    bs2 :: IndexSet

    rng :: MersenneTwister
    function QGomeaMixer(
        f :: Function,
        n :: Int64,
        population :: Vector{QGomeaSolution},
        D :: Matrix{Float64},
        fos :: Vector{Vector{Int64}},
        forced_improvement :: Symbol,
        fos_type :: Symbol,
        crf :: ClusteringReductionFormula,
        repair :: Symbol,
        rng :: MersenneTwister)
        #
        new(f, n, population, D, fos,
            forced_improvement, fos_type, crf, repair,
            Ref(maximum(population)),
            Ref(0), Ref(0), Ref(false),
            collect(1:n), collect(1:n), collect(1:n), collect(1:n),
            IndexSet(n), IndexSet(n),
            rng)
    end

    function QGomeaMixer(
        f :: Function,
        n :: Int64,
        population :: Vector{QGomeaSolution},
        D :: Matrix{Float64},
        fos :: Vector{Vector{Int64}},
        forced_improvement :: Symbol,
        fos_type :: Symbol,
        crf :: ClusteringReductionFormula,
        repair :: Symbol,
        best :: Ref{QGomeaSolution},
        rng :: MersenneTwister)
        #
        best[] = max(best[], maximum(population))
        new(f, n, population, D, fos,
            forced_improvement, fos_type, crf, repair,
            best,
            Ref(0), Ref(0), Ref(false),
            collect(1:n), collect(1:n), collect(1:n), collect(1:n),
            IndexSet(n), IndexSet(n),
            rng)
    end
end

# function mix!(dst :: Vector{T}, donor :: Vector{T}, viewMask, storage :: Vector{Int64}) where {T}
#     dstview = view(dst, viewMask)
#     sort!(dstview)
#     resize!(storage, length(dstview))
#     Base.invpermute!!(dstview, sortperm!(storage, view(donor, viewMask)))
#     dst
# end

# Using Built-In BitSets.
# 2x as fast.
function mix_order!(dst :: Vector{Int64}, donor :: Vector{Int64}, viewMask :: Vector{Int64},
        ds1 :: BitSet, ds2 :: BitSet, storage :: Vector{Int64})
    #
    empty!(ds1); empty!(ds2)
    @inbounds for o in viewMask
        push!(ds1, dst[o])
        push!(ds2, donor[o])
    end
    @inbounds for (nds1,nds2) in zip(ds1, ds2)
        storage[nds2] = nds1
    end
    @inbounds for o in viewMask
        dst[o] = storage[donor[o]]
    end
end

# Using own IndexSet type, somewhat faster due to less checks.
# ~ 4x as fast.
function mix_order!(dst :: Vector{Int64}, donor :: Vector{Int64}, viewMask :: Vector{Int64},
        ds1 :: IndexSet{T}, ds2 :: IndexSet{T}, storage :: Vector{Int64}) where {T}
    #
    empty!(ds1); empty!(ds2)
    @inbounds for o in viewMask
        push!(ds1, convert(UInt64, dst[o]))
        push!(ds2, convert(UInt64, donor[o]))
    end
    for (nds1,nds2) in zip(ds1, ds2)
        storage[nds2] = nds1
    end
    @inbounds for o in viewMask
        dst[o] = storage[donor[o]]
    end
end

function evaluate(f :: Function, perm :: Vector{Int64})
    fitness = f(perm)
    return QGomeaSolution(perm, fitness)
end

function calcD!(pm :: QGomeaMixer)
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

function calcD_original!(pm :: QGomeaMixer)
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

function calcD_original_invd2!(pm :: QGomeaMixer)
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
            δ₂ = 1.0 - entropy(c_ab / pm.n)
            pm.D[i, j] = δ₁ * δ₂
            pm.D[j, i] = pm.D[i, j]
        end
    end
    pm.D
end

function calcD_original_regularized!(pm :: QGomeaMixer)
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

function calcD_random!(pm :: QGomeaMixer, rng :: MersenneTwister)
    for i in 1:pm.n
        for j in i+1:pm.n
            pm.D[i, j] = rand(rng)
            pm.D[j, i] = pm.D[i, j]
        end
    end
    pm.D
end


function findchildrenpairs(parents :: Vector{Int64})
    idx = zero(parents)
    result = Vector{Pair{Int64, Tuple{Int64, Int64}}}()
    for (me, parent) in enumerate(parents)
        if parent == 0
            # Final node, does not have a parent.
            continue
        end
        if idx[parent] == 0
            idx[parent] = me
        else
            push!(result, parent => (idx[parent], me))
        end
    end
    return result
end

function lsmixing(pm :: QGomeaMixer,
    parentChild :: Vector{Pair{Int64, Tuple{Int64, Int64}}};
    shuffle_parent_child = true)
    #
    improved = false
    if shuffle_parent_child
        shuffle!(pm.rng, parentChild)
    end
    for individual in pm.population
        if rand(pm.rng) < 0.9
            continue
        end
        copyto!(pm.mixing_backup, individual.perm)
        for parent_children in parentChild
            parent = parent_children.first
            (child_a, child_b) = parent_children.second
            reverse!(view(pm.mixing_backup, pm.fos[parent]))
            reverse!(view(pm.mixing_backup, pm.fos[child_a]))
            reverse!(view(pm.mixing_backup, pm.fos[child_b]))
            new_fitness = pm.f(pm.mixing_backup)
            if new_fitness > individual.fitness
                copyto!(individual.perm, pm.mixing_backup)
                improved = true
                individual.fitness = new_fitness
            # elseif new_fitness == individual.fitness
            #     copyto!(individual.perm, pm.mixing_backup)
            elseif new_fitness < individual.fitness
                copyto!(pm.mixing_backup, individual.perm)
            end
        end
    end
    return improved
end

function edamixing(sol :: QGomeaSolution, pm :: QGomeaMixer; shuffle_fos=true, donor_fixed :: Union{Nothing, Ref{QGomeaSolution}} = nothing)
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
        #mix!(pm.mixing_backup, donor, s, pm.mixing_perm)
        # If subsets are equal in value, no mixing is going to make
        # a difference
        if all(pm.mixing_backup[i] == donor[i] for i in s)
            continue
        end
        for op in 1:2
            if op === 1
                create_differential_donor!(pm.mixing_virtual_donor, pm.mixing_backup, donor, s, pm.rng,
                    offset=0) #Can also be :rand_difference or :rand_neg1to1
                if pm.repair == :ox
                    mix_position_ox!(pm.mixing_backup, pm.mixing_virtual_donor, s, pm.mixing_perm, pm.mixing_perm_2)
                elseif pm.repair == :pmx
                    mix_position_pmx!(pm.mixing_backup, pm.mixing_virtual_donor, s, pm.mixing_perm, pm.mixing_perm_2)
                else
                    error("Unknown permutation repair approach.")
                end
            elseif op === 2 #
                # This operator is a no-op on |s| == 1, no mixing and evaluation required.
                if length(s) == 1# && length(s) >= round(pm.n/12*10)
                    continue
                end
                mix_order!(pm.mixing_backup, donor, s, pm.bs1, pm.bs2, pm.mixing_perm)
            elseif op === 3
                mix_position_ox!(pm.mixing_backup, donor, s, pm.mixing_perm, pm.mixing_perm_2)
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
    end
    return best_improved, solution_changed
end

function step!(pm :: QGomeaMixer)
    #
    #println([a.fitness for a in sort(pm.population, rev=true)[1:10]])
    # Calculate D depending on the population.
    if pm.fos_type == :distance
        calcD!(pm)
    elseif pm.fos_type == :original
        calcD_original!(pm)
    elseif pm.fos_type == :random
        calcD_random!(pm)
    else
        error("Unknown FOS type $(pm.fos_type)")
    end
    # calcD_original_regularized!(pm)
    # calcD_original_invd2!(pm)
    # Compute FOS
    empty!(pm.fos)
    parent_idx = zeros(Int64, 2*pm.n-1)
    fos_indexset = LCP(pm.D, pm.crf; parent_idx=parent_idx)
    append!(pm.fos, collect(a) for (i,a) in enumerate(fos_indexset))

    improved_any = false
    # Perform LS Mixing.
    # parentChild = findchildrenpairs(parent_idx)
    # improved_any |= lsmixing(pm, parentChild)

    # - Distance based pruning
    # merge_distance = zeros(Float64, 2*pm.n-1)
    # fos_indexset = LCP(pm.D, crf; merge_distance=merge_distance, parent_idx=parent_idx)
    # mdq95 = quantile(merge_distance, 0.999)
    # append!(pm.fos, collect(a) for (i,a) in enumerate(fos_indexset) if merge_distance[i] <= mdq95)

    # append!(pm.fos, [i, j] for i in 1:pm.n for j in (i+1):pm.n)

    # Inv distance
    # fos_indexset = LCP(maximum(pm.D) .- pm.D, crf)
    # append!(pm.fos, collect(a) for a in fos_indexset)

    # Filter it.
    # sort!(pm.fos, by=f->length(f))
    # filter!(f -> (length(f) > 1 && length(f) < round(pm.n/12*10)), pm.fos)
    filter!(f -> (length(f) < round(pm.n/12*10)), pm.fos)
    # filter!(f -> (length(f) > 1), pm.fos)
    # filter!(f -> (length(f) > 1 && length(f) < pm.n), pm.fos)
    # filter!(f -> (length(f) < pm.n), pm.fos)
    # Mixymix.
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

function generate_new_qgomeasolution_random(f :: Function, n :: Int64, rng :: MersenneTwister)
    perm = shuffle!(rng, collect(1:n))
    QGomeaSolution(perm, f(perm))
end

function generate_new_qgomeasolution_bounds(bounds :: Vector{Tuple{Float64, Float64}})
    function generate_new_solution(f :: Function, n :: Int64, rng :: MersenneTwister)
        @assert length(bounds) == n
        perm = invperm(sortperm([rand(rng) * (a[2]-a[1]) + a[1] for a in bounds]))
        QGomeaSolution(perm, f(perm))
    end
    return generate_new_solution
end

function generate_new_qgomeasolution_bounds(bounds :: Tuple{Vector{Float64}, Vector{Float64}})
    function generate_new_solution(f :: Function, n :: Int64, rng :: MersenneTwister)
        @assert length(bounds) == n
        perm = invperm(sortperm([rand(rng) * (ub-lb) + lb for (lb, ub) in zip(bounds[1], bounds[2])]))
        QGomeaSolution(perm, f(perm))
    end
    return generate_new_solution
end

function create_qgomea_mixer(f :: Function, n :: Int64, population_size :: Int64, 
    forced_improvement :: Symbol,
    fos_type :: Symbol,
    crf :: ClusteringReductionFormula,
    permutation_repair :: Symbol,
    rng :: MersenneTwister,
    best :: Union{Nothing, Ref{QGomeaSolution}} = nothing;
    initial_solution_generator :: Function) :: QGomeaMixer
    # Generate initial population.
    population = [initial_solution_generator(f, n, rng) for _ in 1:population_size]
    # Create matrices and FoS vector.
    D = zeros(Float64, (n, n))
    fos = Vector{Vector{Int64}}()
    sizehint!(fos, 2*n)
    if best === nothing
        return QGomeaMixer(f, n, population, D, fos, forced_improvement, fos_type, crf, permutation_repair, rng)
    else
        return QGomeaMixer(f, n, population, D, fos, forced_improvement, fos_type, crf, permutation_repair, best, rng)
    end
end

function optimize_qgomea(rf :: Function, n :: Int64, t=10.0, e=typemax(Int64);
    initial_solution_generator :: Function = generate_new_qgomeasolution_random,
    population_size_base=4, 
    crf=UPGMA(), 
    forced_improvement :: Symbol = :extended,
    fos_type :: Symbol = :distance,
    permutation_repair = :ox, 
    target_fitness :: Union{Nothing, Float64} = nothing)
    #
    time_start = time()
    n_evals = 0

    rng = MersenneTwister()

    fx = wrap_assignment_to_permutation(rf)
    function f(sol :: Vector{Int64})
        n_evals += 1
        return fx(sol)
    end

    next_population_size = population_size_base*2

    initial_mixer = create_qgomea_mixer(f, n, population_size_base, forced_improvement, fos_type, crf, permutation_repair, rng, initial_solution_generator=initial_solution_generator)
    mixers = QGomeaMixer[initial_mixer]
    steps = 0
    last_steps = 0
    best = initial_mixer.best

    while (time() - time_start < t) && (n_evals <= e) && best[].fitness != target_fitness

        for i_mixer in 1:length(mixers)
            # Other steps!
            if mod(steps, 4^(i_mixer-1)) != 0
                break
            end
            if i_mixer == length(mixers)
                last_steps += 1
            end

            step!(mixers[i_mixer])
        end

        # Add new metaheuristic upon having stepped the last one 4 times.
        # Or if all have converged. Oh well.
        if last_steps == 4 || length(mixers) == 0
            last_steps = 0
            push!(mixers, create_qgomea_mixer(f, n, next_population_size, forced_improvement, fos_type, crf, permutation_repair, rng, best, initial_solution_generator=initial_solution_generator))
            next_population_size *= 2
        end
        filter!(f -> !f.converged[], mixers)

        steps += 1
    end

    return (best[].fitness, invperm(best[].perm))
end
