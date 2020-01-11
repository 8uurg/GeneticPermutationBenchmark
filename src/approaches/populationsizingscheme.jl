
"""
A generator that generates the highest index to run for a particular step iteration. 
"""
mutable struct SubpopulationStepcountGenerator
    stack :: Vector{Int64}
    factor :: Int64

    function SubpopulationStepcountGenerator(factor :: Int64)
        return new([1], factor)
    end
end

function step!(g :: SubpopulationStepcountGenerator)
    i = 1
    while i <= length(g.stack) && g.stack[i] == g.factor
        g.stack[i] = 1
        i += 1
    end
    if i > length(g.stack)
        push!(g.stack, 1)
    end
    g.stack[i] += 1
    return i
end

# Add iteration protocol.
Base.iterate(g :: SubpopulationStepcountGenerator) = step!(g), nothing
Base.iterate(g :: SubpopulationStepcountGenerator, :: Nothing) = step!(g), nothing
# Define size.
Base.IteratorSize(g :: SubpopulationStepcountGenerator) = Base.IsInfinite()