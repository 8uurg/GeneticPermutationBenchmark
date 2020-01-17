module GomeaWrapper

import Libdl

# case  0: return( (char *) "Univariate" );
# case  1: return( (char *) "Linkage Tree" );
# case  2: return( (char *) "Multiscale Linkage Neighbors" );
# case  3: return( (char *) "Linkage Trees and Neighbors" );
# case  4: return( (char *) "Filtered Linkage Tree" );
# case  5: return( (char *) "Filtered Multiscale Linkage Neighbors" );
# case  6: return( (char *) "Filtered Linkage Trees and Neighbors" );
FOSStructures = Dict(
    :univariate => 0,
    :linkage_tree => 1,
    :multiscale_linkage_neighbors => 2,
    :linkage_trees_and_neighbors => 3,
    :filtered_linkage_tree => 4,
    :filtered_multiscale_linkage_neighbors => 5,
    :filtered_linkage_trees_and_neighbors => 6
)

"""
    calledInC(functrunk :: Ptr{Cvoid}, g_n :: Int64, numbers :: Ptr{Cdouble}, objective :: Ref{Cdouble}, constraint :: Ref{Cdouble})

Function that is called in C to provide a normal Julia function call instead.

# Arguments
functrunk :: Ptr{Cvoid}     - The void pointer that is turned back into a Julia function, allows for using a closure rather
                              than being limited to a C function pointer.
g_n :: Int64                - The number of parameters in the solution vector `numbers`
numbers :: Ptr{Cdouble}     - Pointer to a Vector of elements.
objective :: Ref{Cdouble}   - Endpoint to write the calculated objective to.
constraint :: Ref{Cdouble} - constraint value - eg. the slack of constraint violated.

"""
function calledInC(functrunk :: Ptr{Cvoid}, g_n :: Int64, numbers :: Ptr{Cdouble}, objective :: Ref{Cdouble}, constraint :: Ref{Cdouble})
    jnumbers = unsafe_wrap(Array, numbers, g_n)
    optimfunc = unsafe_pointer_to_objref(functrunk)::Function
    obj, feas = optimfunc(jnumbers)
    unsafe_store!(objective, obj)
    unsafe_store!(constraint, feas)
    nothing
end

"""
    optimizeGOMEA(func :: Function, n :: Int64; maxTime :: Int64 = 100, maxEvaluations = -1)

Optimize a function using GOMEA

# Arguments
func :: Function        - A function that accepts a vector of `n` Float64 values and returns a tuple (objective, constraint)
n :: Int64              - The amount of parameters that needs to be optimized over
maxTime :: Int64        - The maximum amount of time optimization may take (in milliseconds), -1 for unlimited.
maxEvaluations :: Int64 - The maximum amount of evaluations of func, -1 for unlimited.
fosStructure :: Symbol  - One of :univariate, :linkage_tree, :multiscale_linkage_neighbors, :linkage_trees_and_neighbors,
                          :filtered_linkage_tree, :filtered_multiscale_linkage_neighbors, :filtered_linkage_trees_and_neighbors
rkeys_until :: Int64    - The index up to which a random key encoding is used.
reencode :: Bool        - Whether reencoding of the area on which random key encoding is enabled.
!!! danger Potential incorrect behavior
    Reencoding will not reevaluate the solution. If you are not using a random key
    encoding, enabling this will cause inconsistencies in the final solution.
rescale :: Bool         - Whether randomly rescaling is enabled.

# Returns
A tuple (objective, constraint, elite solution)

"""
function optimizeGOMEA(func :: Function, n :: Int64; 
    maxTime :: Int64 = 100, maxEvaluations :: Int64 = -1, # Execution limits.
    fosStructure :: Symbol = :univariate, 
    rkeys_until :: Int64 = n,
    bounds :: Union{Tuple{Float64, Float64}, 
                    Tuple{Vector{Float64}, Vector{Float64}},
                    Vector{Tuple{Float64, Float64}}} = (0.0, 1.0), 
    reencode :: Bool = false, rescale :: Bool = true, translate :: Bool = false,
    performs_ls :: Bool = false,
    dependency_matrix_approach :: Symbol = :original )
    # Safety checks, assuming 1-indexing
    @assert rkeys_until <= n 

    # Declare some memory for return values and elite solution.
    objective = Ref{Cdouble}(0)
    constraint = Ref{Cdouble}(0)
    best_solution = zeros(Cdouble, n)

    fos = FOSStructures[fosStructure]
    eva = maxEvaluations
    mil = maxTime

    lowerbounds :: Vector{Float64} = [0.0]
    upperbounds :: Vector{Float64} = [1.0]
    boundsmode :: Int8 = 0

    dependency_matrix_approach_char = 0
    if dependency_matrix_approach == :avg_distance
        dependency_matrix_approach_char = 1
    elseif dependency_matrix_approach == :random
        dependency_matrix_approach_char = 2
    elseif dependency_matrix_approach == :inv_avg_distance
        dependency_matrix_approach_char = 3
    end

    if typeof(bounds) == Tuple{Float64, Float64}
        boundsmode = 0
        lowerbounds = [0.0]
        upperbounds = [1.0]
    elseif typeof(bounds) == Tuple{Vector{Float64}, Vector{Float64}}
        # Bounds bounds checks.
        @assert length(bounds[1]) == n
        @assert length(bounds[2]) == n
        boundsmode = 1
        lowerbounds = bounds[1]
        upperbounds = bounds[2]
    elseif typeof(bounds) == Vector{Tuple{Float64, Float64}}
        @assert length(bounds) == n
        boundsmode = 1
        lowerbounds = [bound[1] for bound in bounds]
        upperbounds = [bound[2] for bound in bounds]
    else
        error("This should never happen by typecheck.")
    end

    __gomea = Libdl.dlopen(joinpath(@__DIR__, "..", "deps", "GOMEA"))
    __optim = Libdl.dlsym(__gomea, "optimizeGOMEA")

    ccalled = @cfunction(calledInC, Cvoid, (Ptr{Cvoid}, Int64, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}))
    ccall(__optim, Cvoid,
       (Any, Ptr{Cvoid}, Int64, Int64, Int64, Int64,       Int64,     Bool,    Bool,      Bool,     Cchar, Ptr{Cdouble}, Ptr{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Vector{Cdouble}},       Cchar,                           Cchar),
        func,   ccalled,     n,   fos,   eva,   mil, rkeys_until, reencode, rescale, translate,boundsmode,  lowerbounds,  upperbounds,    objective,   constraint,        best_solution, performs_ls, dependency_matrix_approach_char)

    (convert(Float64, objective[]), convert(Float64, constraint[]), convert(Vector{Float64}, best_solution))
end

"""
    objectiveOnly(f :: Function)

Turn a function that only returns an objective value (and not also a constraint value)
into one that is accepted by optimizeGOMEA.

# Arguments
f :: Function - A function that accepts a vector of doubles.

# Returns
A function that returns both the objective value from f, and 0.0 for the constraint value.
"""
function objectiveOnly(f :: Function)
    function func(vec :: Vector{Float64})
        f(vec), 0.0
    end
end

function sortedness(params :: Vector{Float64})
    obj = sum(Float64, sortperm(params, rev=true) .* (1:length(params)))
    -obj, 0.0
end

function sortedness_pairs(params :: Vector{Float64})
    perm = sortperm(params)
    obj = sum(i+1==j for (i,j) in zip(perm, Iterators.drop(perm, 1)))
    obj, 0.0
end

end # module
