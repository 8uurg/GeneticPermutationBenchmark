# This file contains some benchmark permutation functions related to 'sortedness'
# !!!note All of the functions are implemented as **maximization** problems.

"""
    sorted_inversion(perm :: Vector{Int64})

Count the number of uninverted pairs in the given permutation `perm`.
Should be maximized.

For example in the case of `[3, 1, 2]` contains three pairs:
    [(3, 1), (1, 2), (3, 2)]
        >       <       >
From these pairs one is in the ascending order (<), and therefore the corresponding
count is 1.

```jldoctest
julia> sorted_inversion([3, 1, 2])
1
```

Optimum is the sorted sequence, [1, ..., n],
as can be obtained by `collect(n:-1:1)` for a given `n`.
"""
function sorted_inversion(perm :: Vector{Int64})
    score = 0
    @inbounds for i in 1:length(perm)
        for j in (i+1):length(perm)
            score += ifelse(perm[i] < perm[j], 1, 0)
        end
    end
    return score
end

"""
    sorted_inversion_rev(perm :: Vector{Int64})

Count the number of inverted pairs in the given permutation `perm`.
Should be maximized.

For example in the case of `[3, 1, 2]` contains three pairs:
    [(3, 1), (1, 2), (3, 2)]
        >       <       >
From these pairs two are in the reverse order (>), and therefore the corresponding
count is 2.

```jldoctest
julia> sorted_inversion_rev([3, 1, 2])
2
```

Optimum is reverse of the sorted sequence, [n, ..., 1],
as can be obtained by `collect(n:-1:1)` for a given `n`.
"""
function sorted_inversion_rev(perm :: Vector{Int64})
    score = 0
    @inbounds for i in 1:length(perm)
        for j in (i+1):length(perm)
            score += ifelse(perm[i] > perm[j], 1, 0)
        end
    end
    return score
end

"""
    sorted_sequential_pairs(perm :: Vector{Int64})

Count the number of sequential elements (a, b) in the given permutation `perm`
where the second (b) is the successor to the first (a). Eg. count the cases in which πᵢ₊₁ = πᵢ + 1 
Should be maximized.

For example in the case of `[3, 1, 2]` this does not hold for (3, 1), but does hold for (1, 2).
Therefore the count is 1.

```jldoctest
julia> sorted_sequential_pairs([3, 1, 2])
1
```

Optimum is the sorted sequence, [1, ..., n],
as can be obtained by `collect(n:-1:1)` for a given `n`.
"""
function sorted_sequential_pairs(perm :: Vector{Int64})
    score = 0
    @inbounds for i in 2:length(perm)
        score += ifelse(perm[i] - 1 == perm[i-1], 1, 0)
    end
    return score
end

"""
    sorted_sequential_pairs_rev(perm :: Vector{Int64})

Count the number of sequential elements (a, b) in the given permutation `perm`
where the first (a) is the successor to the second (b). Eg. count the cases in which πᵢ = πᵢ₊₁ + 1 
Should be maximized.

For example in the case of `[1, 3, 2]` this does not hold for (1, 3), but does hold for (3, 2).
Therefore the count is 1.

```jldoctest
julia> sorted_sequential_pairs([1, 3, 2])
1
```
    
Optimum is reverse of the sorted sequence, [n, ..., 1],
as can be obtained by `collect(n:-1:1)` for a given `n`.
"""
function sorted_sequential_pairs_rev(perm :: Vector{Int64})
    score = 0
    @inbounds for i in 2:length(perm)
        score += ifelse(perm[i] + 1 == perm[i-1], 1, 0)
    end
    return score
end

"""
    sorted_sequential_inversion(perm :: Vector{Int64})

Count the number of sequential elements (a, b) in the given permutation `perm`
where the second (b) is greater than the first (a). Eg. count the cases in which πᵢ₊₁ > πᵢ 
Should be maximized.

For example in the case of `[3, 1, 2]` this does not hold for (3, 1), but does hold for (1, 2).
Therefore the count is 1.

```jldoctest
julia> sorted_sequential_inversion([3, 1, 2])
1
```

Optimum is the sorted sequence, [1, ..., n],
as can be obtained by `collect(n:-1:1)` for a given `n`.
"""
function sorted_sequential_inversion(perm :: Vector{Int64})
    score = 0
    for i in 1:(length(perm)-1)
        score += ifelse(perm[i] < perm[i+1], 1, 0)
    end
    return score
end

"""
    sorted_sequential_inversion_rev(perm :: Vector{Int64})

Count the number of sequential elements (a, b) in the given permutation `perm`
where the first (a) is greater than the second (b). Eg. count the cases in which πᵢ > πᵢ₊₁ 
Should be maximized.

For example in the case of `[1, 3, 2]` this does not hold for (1, 3), but does hold for (3, 2).
Therefore the count is 1.

```jldoctest
julia> sorted_sequential_inversion_rev([1, 3, 2])
1
```
    
Optimum is reverse of the sorted sequence, [n, ..., 1],
as can be obtained by `collect(n:-1:1)` for a given `n`.
"""
function sorted_sequential_inversion_rev(perm :: Vector{Int64})
    score = 0
    for i in 1:(length(perm)-1)
        score += ifelse(perm[i] > perm[i+1], 1, 0)
    end
    return score
end

struct BIVInstance{F <: Function}
    func :: F
    n :: Int64
end

function bb_wrap_biv(instance :: BIVInstance)
    function evaluate(permutation :: Vector{Int64}) :: Float64
        return instance.func(permutation)
    end
    return evaluate
end

#
function generate_good_fos_sorted_inversion(n :: Int64)
    return vcat([[i] for i in 1:n], [[i,i+1] for i in 1:n-1])
end