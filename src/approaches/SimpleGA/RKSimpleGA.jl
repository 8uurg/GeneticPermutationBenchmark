using Random
include("../populationsizingscheme.jl")

function wrap_rkeys_to_permutation(f :: Function)
    prm = collect(1:100)
    function ev(assignment :: Vector{Float64})
        resize!(prm, length(assignment))
        sortperm!(prm, assignment)
        return f(prm)
    end
    return ev
end

mutable struct RKSimpleGASolution
    perm :: Vector{Float64}
    fitness :: Float64
end
Base.isless(a :: RKSimpleGASolution, b :: RKSimpleGASolution) = isless(a.fitness, b.fitness)
Base.isequal(a :: RKSimpleGASolution, b :: RKSimpleGASolution) = a.perm == b.perm

struct RKSimpleGA
    f :: Function
    n :: Int64
    population :: Vector{RKSimpleGASolution}
    # Stats & info
    best :: Ref{RKSimpleGASolution}

    generations :: Ref{Int64}
    generations_no_improvement :: Ref{Int64}

    converged :: Ref{Bool}

    # Placeholders
    offspring1 :: RKSimpleGASolution
    offspring2 :: RKSimpleGASolution

    rng :: MersenneTwister

    function RKSimpleGA(f :: Function, n :: Int64, population :: Vector{RKSimpleGASolution},
            rng :: MersenneTwister)
        new(f, n, population,
            # Stats & Info
            Ref(maximum(population)), 
            Ref(0), Ref(0), Ref(false),
            # Placeholders
            RKSimpleGASolution(collect(1:n), typemin(Float64)),
            RKSimpleGASolution(collect(1:n), typemin(Float64)),
            rng)
    end

    function RKSimpleGA(f :: Function, n :: Int64, population :: Vector{RKSimpleGASolution},
        best :: Ref{RKSimpleGASolution}, rng :: MersenneTwister)
        best[] = max(best[], maximum(population))
        new(f, n, population,
            best, Ref(0), Ref(0), Ref(false),
            # Placeholders
            RKSimpleGASolution(collect(1:n), typemin(Float64)),
            RKSimpleGASolution(collect(1:n), typemin(Float64)),
            rng)
    end
end

function uniform_rk_crossover!(offspring1 :: Vector{Float64}, offspring2 :: Vector{Float64},
    parent1 :: Vector{Float64}, parent2 :: Vector{Float64}, rng :: MersenneTwister, p :: Float64 = 0.5)
    #
    n = length(parent1)
    # Copy over
    copy!(offspring1, parent1)
    copy!(offspring2, parent2)
    # Crossover
    @inbounds for i in 1:n
        if rand(rng) < p
            offspring1[i], offspring2[i] = offspring2[i], offspring1[i]
        end
    end
end
function compare_swap!(a :: RKSimpleGASolution, b :: RKSimpleGASolution)
    # Note: <= can be used as well, but you might want to return false instead
    #       when the fitness is equal.
    if a.fitness < b.fitness
        a.perm, b.perm = b.perm, a.perm
        a.fitness, b.fitness = b.fitness, a.fitness
        return true
    end
    return false
end

function elisist_crossover!(ga :: RKSimpleGA, dst1 :: RKSimpleGASolution, dst2 :: RKSimpleGASolution)
    improved_a_solution = false
    improved_best = false
    # Create offspring
    uniform_rk_crossover!(ga.offspring1.perm, ga.offspring2.perm, dst1.perm, dst2.perm, ga.rng, 0.5)
    # Evaluate offspring
    ga.offspring1.fitness = ga.f(ga.offspring1.perm)
    ga.offspring2.fitness = ga.f(ga.offspring2.perm)
    # Lazy 'sorting' net, best two end up in dst1 and dst2, remainder ends up
    # assigned to offspring, to be reused as memory for another crossover.
    improved_a_solution |= compare_swap!(dst1, ga.offspring1)
    improved_a_solution |= compare_swap!(dst1, ga.offspring2)
    improved_a_solution |= compare_swap!(dst2, ga.offspring1)
    improved_a_solution |= compare_swap!(dst2, ga.offspring2)
    # Update best
    if dst1.fitness > ga.best[].fitness
        ga.best[] = dst1
        improved_best = true
    end
    if dst2.fitness > ga.best[].fitness
        ga.best[] = dst2
        improved_best = true
    end
    # 
    return improved_a_solution, improved_best
end

function step!(ga :: RKSimpleGA)
    # 
    improved_a_solution = false
    improved_best = false
    # Shuffle the population.
    shuffle!(ga.population)
    # Perform crossovers
    for i in 1:2:length(ga.population)-1
        improved_a_solution_cc, improved_best_cc =
            elisist_crossover!(ga, ga.population[i], ga.population[i+1])
        improved_a_solution |= improved_a_solution_cc
        improved_best |= improved_best_cc
    end
    # Update statistics
    ga.generations[] += 1

    if !improved_a_solution
        ga.generations_no_improvement[] += 1

        NIS_convergence_threshold = 10 + floor(Int64, 20 * log10(length(ga.population)))
        if ga.generations_no_improvement[] > NIS_convergence_threshold
            ga.converged[] = true
        end
    else
        ga.generations_no_improvement[] = 0
    end
end

function generate_new_RKsimplegasolution_random(f :: Function, n :: Int64, rng :: MersenneTwister)
    perm = rand(rng, Float64, n)
    RKSimpleGASolution(perm, f(perm))
end

function create_RKsimplega(f :: Function, n :: Int64, population_size :: Int64,
    rng :: MersenneTwister,
    best :: Union{Nothing, Ref{RKSimpleGASolution}} = nothing;
    initial_solution_generator :: Function) :: RKSimpleGA
    # Generate initial population.
    population = [initial_solution_generator(f, n, rng) for _ in 1:population_size]
    #
    if best === nothing
        return RKSimpleGA(f, n, population, rng)
    else
        return RKSimpleGA(f, n, population, best, rng)
    end
end

function optimize_rksimplega(rf :: Function, n :: Int64, t=10.0, e=typemax(Int64);
    initial_solution_generator :: Function = generate_new_RKsimplegasolution_random,
    population_size_base=16, population_sizing_step_factor :: Int64 = 8, target_fitness :: Union{Nothing, Float64} = nothing)
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

    initial_mixer = create_RKsimplega(f, n, population_size_base, rng, initial_solution_generator=initial_solution_generator)
    mixers = RKSimpleGA[initial_mixer]
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
            push!(mixers, create_RKsimplega(f, n, next_population_size, rng, best, initial_solution_generator=initial_solution_generator))
            # println("Created new population of size $(next_population_size)")
            step!(mixers[end])
            next_population_size *= 2
        end
        filter!(f -> !f.converged[], mixers)

        steps += 1
    end

    return (best[].fitness, sortperm(best[].perm))
end
