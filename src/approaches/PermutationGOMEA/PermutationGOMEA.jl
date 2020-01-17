using GomeaWrapper

function wrap_objective(f :: Function)
    function to_multiple(x :: Vector{Float64})
        return f(x), 0.0
    end
    return to_multiple
end

function optimize_permutationgomea(func :: Function, n :: Int64, t :: Int64, e :: Int64, 
    fos_type :: Symbol, reencode=true, rescale=true)
    #
    fosStructure = :linkage_tree
    if fos_type == :univariate
        fosStructure = :univariate
    end
    optimizeGOMEA(wrap_objective(func), n, t, e, fosStructure, reencode=reencode, rescale=rescale)
end