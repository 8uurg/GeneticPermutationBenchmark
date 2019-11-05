using Random
include("./recombinationops.jl")

mutable struct IPSimpleGASolution
    perm :: Vector{Int64}
    fitness :: Float64
end

struct IPSimpleGA{O <: CrossoverOperator}
    f :: Function,
    n :: Int64,
    population :: Vector{IPSimpleGASolution}
    crossover :: O
    # Stats & info
    best :: Ref{IPSimpleGASolution}

    generations :: Ref{Int64}
    generations_no_improvement :: Ref{Int64}

    # Placeholders
    offspring1 :: IPSimpleGASolution
    offspring2 :: IPSimpleGASolution

    function IPSimpleGA{O <: CrossoverOperator}(f :: Function, n :: Int64, population :: Vector{IPSimpleGASolution},
        crossover :: O)
        new(f, n, population, crossover, 
            # Stats & Info
            Ref(maximum(population)), 
            Ref(0), Ref(0), 
            # Placeholders
            IPSimpleGASolution(collect(1:n), typemin(Float64)),
            IPSimpleGASolution(collect(1:n), typemin(Float64)))
    end

    function IPSimpleGA{O <: CrossoverOperator}(f :: Function, n :: Int64, population :: Vector{IPSimpleGASolution},
        crossover :: O, best :: Ref{IPSimpleGASolution})
        best[] = max(best[], maximum(population))
        new(f, n, population, crossover, 
            best, Ref(0), Ref(0),
            # Placeholders
            IPSimpleGASolution(collect(1:n), typemin(Float64)),
            IPSimpleGASolution(collect(1:n), typemin(Float64)))
    end
end

function compare_swap!(a :: IPSimpleGASolution, b :: IPSimpleGASolution)
    # Note: <= can be used as well, but you might want to return false instead
    #       when the fitness is equal.
    if a.fitness < b.fitness
        a.perm, b.perm = b.perm, a.perm
        a.fitness, b.fitness = b.fitness, a.fitness
        return true
    end
    return false
end

function elisist_crossover!(ga :: IPSimpleGA{O}, dst1 :: IPSimpleGASolution, dst2 :: IPSimpleGASolution) where {O}
    improved_a_solution = false
    improved_best = false
    # Create offspring
    crossover!(ga.crossover, ga.offspring1.perm, ga.offspring2.perm, dst1.perm, dst2.perm)
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

function step!(ga :: IPSimpleGA)
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
        ga.generations_no_improvement += 1
    end
end

function generate_new_ipsimplegasolution_random(f :: Function, n :: Int64)
    perm = shuffle!(collect(1:n))
    IPSimpleGASolution(perm, f(perm))
end

function create_ipsimplega(f :: Function, n :: Int64, population_size :: Int64, crossover :: O,
    best :: Union{Nothing, Ref{IPSimpleGASolution}} = nothing;
    initial_solution_generator :: Function) :: IPSimpleGA where {O <: CrossoverOperator}
    # Generate initial population.
    population = [initial_solution_generator(f, n) for _ in 1:population_size]
    #
    if best === nothing
        return IPSimpleGA(f, n, population, crossover)
    else
        return IPSimpleGA(f, n, population, crossover, best)
    end
end

function optimize_ipsimplega(crossover :: O, f :: Function, n :: Int64, t=10.0;
    initial_solution_generator :: Function = generate_new_ipsimplegasolution_random,
    population_size_base=4, target_fitness :: Union{Nothing, Float64} = nothing) where {O <: CrossoverOperator}
    #
    time_start = time()

    next_population_size = population_size_base*2

    initial_mixer = create_ipsimplega(f, n, population_size_base, crossover, initial_solution_generator=initial_solution_generator)
    mixers = IPSimpleGA[initial_mixer]
    steps = 0
    last_steps = 0
    best = initial_mixer.best

    while (time() - time_start < t) && best[].fitness != target_fitness

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
            push!(mixers, create_ipsimplega(crossover, f, n, next_population_size, best, initial_solution_generator=initial_solution_generator))
            next_population_size *= 2
        end
        filter!(f -> !f.converged[], mixers)

        steps += 1
    end

    return (best[].fitness, best[].perm)
end
