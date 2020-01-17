using GomeaWrapper

function wrap_objective(f :: Function, n :: Int64)
    perm = collect(1:n)
    function to_multiple(x :: Vector{Float64})
        sortperm!(perm, x)
        return f(perm), 0.0
    end
    return to_multiple
end

function optimize_permutationgomea(func :: Function, n :: Int64, t :: Float64, e :: Int64 = -1; 
    fos_type :: Symbol = :original, reencode=true, rescale=true)
    #
    fosStructure = :linkage_tree
    if fos_type == :univariate
        fosStructure = :univariate
    end
    conv_t = -1
    if t != Inf
        convert(Int64, t*1000)
    end
    obj, constr, best = optimizeGOMEA(wrap_objective(func, n), n; maxTime=conv_t, maxEvaluations=e, fosStructure=fosStructure, dependency_matrix_approach=fos_type, reencode=reencode, rescale=rescale)
    (obj, sortperm(best))
end