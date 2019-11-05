using Random

mutable struct RKSimpleGASolution
    perm :: Vector{Float64}
    fitness :: Float64
end

struct RKSimpleGA
    f :: Function
    n :: Int64
    population :: Vector{RKSimpleGASolution}
    # Stats & info
    best :: Ref{RKSimpleGASolution}

    generations :: Ref{Int64}
    generations_no_improvement :: Ref{Int64}

    # Placeholders
    offspring1 :: RKSimpleGASolution
    offspring2 :: RKSimpleGASolution

    function RKSimpleGA(f :: Function, n :: Int64, population :: Vector{RKSimpleGASolution})
        new(f, n, population,
            # Stats & Info
            Ref(maximum(population)), 
            Ref(0), Ref(0), 
            # Placeholders
            RKSimpleGASolution(collect(1:n), typemin(Float64)),
            RKSimpleGASolution(collect(1:n), typemin(Float64)))
    end

    function RKSimpleGA(f :: Function, n :: Int64, population :: Vector{RKSimpleGASolution},
        best :: Ref{RKSimpleGASolution})
        best[] = max(best[], maximum(population))
        new(f, n, population,
            best, Ref(0), Ref(0),
            # Placeholders
            RKSimpleGASolution(collect(1:n), typemin(Float64)),
            RKSimpleGASolution(collect(1:n), typemin(Float64)))
    end
end

function uniform_rk_crossover!(offspring1 :: Vector{Int64}, offspring2 :: Vector{Int64},
    parent1 :: Vector{Int64}, parent2 :: Vector{Int64}; p = 0.5)
    #
    n = length(parent1)
    # Copy over
    copy!(offspring1, parent1)
    copy!(offspring2, parent2)
    # Crossover
    @inbounds for i in 1:n
        if rand() < p
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

function elisist_crossover!(ga :: RKSimpleGA{O}, dst1 :: RKSimpleGASolution, dst2 :: RKSimpleGASolution) where {O}
    improved_a_solution = false
    improved_best = false
    # Create offspring
    uniform_rk_crossover!(ga.offspring1.perm, ga.offspring2.perm, dst1.perm, dst2.perm)
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
        ga.generations_no_improvement += 1
    end
end

function generate_new_RKsimplegasolution_random(f :: Function, n :: Int64)
    perm = shuffle!(collect(1:n))
    RKSimpleGASolution(perm, f(perm))
end

function create_RKsimplega(f :: Function, n :: Int64, population_size :: Int64,
    best :: Union{Nothing, Ref{RKSimpleGASolution}} = nothing;
    initial_solution_generator :: Function) :: RKSimpleGA
    # Generate initial population.
    population = [initial_solution_generator(f, n) for _ in 1:population_size]
    #
    if best === nothing
        return RKSimpleGA(f, n, population)
    else
        return RKSimpleGA(f, n, population, best)
    end
end

function optimize_RKsimplega(f :: Function, n :: Int64, t=10.0;
    initial_solution_generator :: Function = generate_new_RKsimplegasolution_random,
    population_size_base=4, target_fitness :: Union{Nothing, Float64} = nothing)
    #
    time_start = time()

    next_population_size = population_size_base*2

    initial_mixer = create_RKsimplega(f, n, population_size_base, initial_solution_generator=initial_solution_generator)
    mixers = RKSimpleGA[initial_mixer]
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
            push!(mixers, create_RKsimplega(f, n, next_population_size, best, initial_solution_generator=initial_solution_generator))
            next_population_size *= 2
        end
        filter!(f -> !f.converged[], mixers)

        steps += 1
    end

    return (best[].fitness, best[].perm)
end
