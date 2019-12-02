using Random

"""
    random_remap(f :: Function, n :: Int64)

Take a black-box function `f` that accepts a permutation as input and create a new black-box
function `rex` that remaps the input permutation using a random mapping.
This remapped solution is then evaluated using `f` and its corresponding objective value is returned.
"""
function random_remap(f :: Function, n :: Int64)
    mapping = shuffle!(collect(1:n))
    placeholder = shuffle!(collect(1:n))
    function rex(perm :: Vector{Int64})
        for i in 1:length(perm)
            placeholder[i] = mapping[perm[i]]
        end
        return f(placeholder)
    end
    return rex
end

"""
    random_remap_and_give_mapping(f :: Function, n :: Int64)

This function does the same thing as random_remap, but in addition the mapping is returned as well.
"""
function random_remap_and_give_mapping(f :: Function, n :: Int64)
    mapping = shuffle!(collect(1:n))
    placeholder = shuffle!(collect(1:n))
    function rex(perm :: Vector{Int64})
        for i in 1:length(perm)
            placeholder[i] = mapping[perm[i]]
        end
        return f(placeholder)
    end
    return rex, mapping
end
