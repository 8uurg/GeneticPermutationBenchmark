include("./IndexSet.jl")
using LinearAlgebra
using Random

abstract type ClusteringReductionFormula end

# Proven to be equivalent to the global strategy.
struct UPGMA <: ClusteringReductionFormula end
struct WPGMA <: ClusteringReductionFormula end
struct SingleLinkage <: ClusteringReductionFormula end
# Not proven, but convex.
struct CompleteLinkage <: ClusteringReductionFormula end

@inline function reduction(::UPGMA, lCₖ, lCᵢ, lCⱼ, DCₖᵢ, DCₖⱼ)
    lCᵢ/(lCᵢ + lCⱼ) * DCₖᵢ + lCⱼ/(lCᵢ + lCⱼ) * DCₖᵢ
end

@inline function reduction(::WPGMA, Cₖ, Cᵢ, Cⱼ, DCₖᵢ, DCₖⱼ)
    (DCₖᵢ + DCₖⱼ) / 2
end

@inline function reduction(::SingleLinkage, Cₖ, Cᵢ, Cⱼ, DCₖᵢ, DCₖⱼ)
    min(DCₖᵢ, DCₖⱼ)
end

@inline function reduction(::CompleteLinkage, Cₖ, Cᵢ, Cⱼ, DCₖᵢ, DCₖⱼ)
    max(DCₖᵢ, DCₖⱼ)
end

function nnmasked(M :: Matrix, j :: Int64, rowmask, rng :: MersenneTwister, randomized :: Bool)
    am = 1
    cx = 0
    cmin = typemax(Float64)
    for (ismasked, (i, v)) in zip(rowmask, enumerate(view(M, :, j)))
        if ismasked && v <= cmin && i != j
            if randomized
                if v < cmin
                    cx = 0
                end
                cx += 1
                if rand(rng, 1:cx) == 1
                    cmin = v
                    am = i
                end
            else
                if v < cmin
                    cmin = v
                    am = i
                end
            end
        end
    end
    am
end

idxunit(x, setsize) = begin
     idxset = IndexSet(setsize);
     push!(idxset, x);
     idxset
 end

function LCP!!(D :: Matrix,
               FoS :: Vector{IndexSet{V}},
               crf :: T,
               rng :: MersenneTwister;
               merge_distance :: Union{Nothing, Vector{Float64}} = nothing,
               parent_idx :: Union{Nothing, Vector{Int64}} = nothing,
               randomized :: Bool) where
               {T <: ClusteringReductionFormula, V}
    #
    @assert issymmetric(D) "D is required to be symmetric."

    n = size(D, 1)
    idxsetsize = _div64(n) + 1
    inuse = trues(n)

    if merge_distance !== nothing
        empty!(merge_distance)
        sizehint!(merge_distance, 2n-1)
        resize!(merge_distance, n)
        fill!(merge_distance, zero(Float64))
    end

    fos_idx = zeros(Int64, n)
    if parent_idx !== nothing
        empty!(parent_idx)
        sizehint!(parent_idx, 2n-1)
        resize!(parent_idx, n)
        fill!(parent_idx, zero(Float64))
        fos_idx .= 1:n
    end
    # Initialize FoS and mpm
    empty!(FoS)
    mpm = Vector{IndexSet{idxsetsize}}(undef, n)
    mpm_size = Vector{Int64}(undef, n)
    sizehint!(FoS, n*2-1)
    resize!(FoS, n)
    for i in 1:n
        mpm[i] = idxunit(i, n)
        mpm_size[i] = 1
        FoS[i] = idxunit(i, n)
    end
    # Set up initial chain.
    NN_chain = Vector{Int64}(undef, 0)
    mpm_s = length(mpm)
    sizehint!(NN_chain, n)
    done = false

    while !done
        # Build chain
        if length(NN_chain) == 0
            # Pick a random element to start off with.
            push!(NN_chain, rand(findall(inuse)))
        end
        while length(NN_chain) < 3
            @inbounds push!(NN_chain, nnmasked(D, NN_chain[end], inuse, rng, randomized))
        end
        # Until reciprocal point found.
        @inbounds while NN_chain[end] != NN_chain[end-2]

            if length(NN_chain) == n
                nxt = NN_chain[end-1]
                push!(NN_chain, nxt)
                break;
            end
            nxt = nnmasked(D, NN_chain[end], inuse, rng, randomized)

            # Make a tiebreaker go in favour of a reciprocal point.
            @inbounds if D[NN_chain[end], nxt] == D[NN_chain[end], NN_chain[end-1]]
                @inbounds nxt = NN_chain[end-1]
            end
            push!(NN_chain, nxt)
        end
        # Merge these two together. r0 will represent the new cluster.
        @inbounds r0, r1 = minmax(NN_chain[end-1], NN_chain[end])
        @inbounds inuse[r0] = false
        @inbounds inuse[r1] = false
        # Shrink the chain
        pop!(NN_chain); pop!(NN_chain); pop!(NN_chain);
        # Update D, lowest item in set has its column and row replaced
        # by that of the new set.
        # view(D, inuse, r0) .= reduction.(Ref(crf),
        #     view(mpm, inuse),
        #     Ref(mpm[r0]), Ref(mpm[r1]),
        #     view(D, r0, inuse), view(D, r1, inuse))
        # view(D, r0, inuse) .= view(D, inuse, r0)
        if merge_distance !== nothing
            push!(merge_distance, D[r0, r1])
        end
        for i in 1:length(mpm)
            @inbounds if inuse[i]
                @inbounds rd = reduction(crf, mpm_size[i], mpm_size[r0], mpm_size[r1], D[r0, i], D[r1, i])
                @inbounds D[i, r0] = rd
                @inbounds D[r0, i] = rd
            end
        end
        mpm_s -= 1
        @inbounds inuse[r0] = true
        # Update set of orders.
        @inbounds append!(mpm[r0], mpm[r1])
        @inbounds mpm_size[r0] += mpm_size[r1]
        #
        push!(FoS, IndexSet(mpm[r0]))
        if parent_idx !== nothing
            parent_idx[fos_idx[r0]] = length(FoS)
            parent_idx[fos_idx[r1]] = length(FoS)
            push!(parent_idx, 0)
            fos_idx[r0] = length(FoS)
        end

        if mpm_s == 1
            done = true
        end
    end
end

function LCP(D :: Matrix,
        crf :: T,
        rng :: MersenneTwister = MersenneTwister();
        merge_distance :: Union{Nothing, Vector{Float64}} = nothing,
        parent_idx :: Union{Nothing, Vector{Int64}} = nothing,
        randomized :: Bool = true) where {T <: ClusteringReductionFormula, V <: Any}
    n = size(D, 1)
    idxsetsize = _div64(n) + 1
    FoS = Vector{IndexSet{idxsetsize}}(undef, 2n-1)
    LCP!!(copy(D), FoS, crf, rng;
        merge_distance=merge_distance,
        parent_idx=parent_idx,
        randomized = randomized)
    return FoS
end

function export_tree_dot(FoS :: Vector{IndexSet{V}}, parent_idx :: Vector{Int64}) where {V}
    """
    digraph {
        $(join([
            string(i, "[label=\"$(join([string(i) for i in FoS[i]], ", "))\"];\n")
        for i in 1:length(FoS)]))
        $(join([
            string(i, " -> ", j, ";\n") 
            for (i,j) in enumerate(parent_idx) if j != 0]))
    }
    """
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

##
function rewrite_by_swap_fos(fos :: Vector{Vector{Int64}},
                             child_pairs :: Vector{Pair{Int64, Tuple{Int64, Int64}}},
                             D :: Matrix{Float64})
    # Constants
    n = size(D,1)
    k = length(fos)
    swap_count = 0
    # Swaps are performed by renaming (eg. lookup table)
    alias = collect(1:n)
    # ΣΣD(A) for A ∈ fos
    DF = zeros(Float64, k)
    # Storage for D(x, A), D(x, B), D(x, A), D(y, B)
    DxF = zeros(Float64, n, k)
    # @inbounds if all(length(fos[i]) == 1 && fos[i][1] == i for i in 1:n)
        # Only works for linkage tree starting with [1], [2], ..., [n-1], [n]
        # DxF[1:n, 1:n] .= D
    # else
    # Works for any Linkage Tree
    @inbounds for (x, f1) in enumerate(fos)
        if length(f1) == 1
            for (y, f2) in enumerate(fos)
                if length(f2) == 1
                    DxF[x, y] = D[first(f1), first(f2)]
                end
            end
        end
    end
    # end
    # Required: Child pairs are ordered such that children occur before parents.
    @inbounds for (fos_idx_a_cup_b, (fos_idx_a, fos_idx_b)) in child_pairs
        # Grab values.
        subset_a = fos[fos_idx_a]
        D_A = DF[fos_idx_a]
        subset_b = fos[fos_idx_b]
        D_B = DF[fos_idx_b]
        # Compute ΣD(a, B)
        for elem_a in subset_a
            DxF[alias[elem_a], fos_idx_b] = sum(D[alias[elem_a], alias[elem_b]] for elem_b in subset_b)
        end
        # Compute ΣD(b, A)
        for elem_b in subset_b
            DxF[alias[elem_b], fos_idx_a] = sum(D[alias[elem_b], alias[elem_a]] for elem_a in subset_a)
        end
        # Compute ΣΣD(A, B)
        # Currently using subset A. Can also use subset B.
        # Maybe pick the smallest of the two. Though I doubt it matters much.
        D_AB = sum(DxF[alias[elem_a], fos_idx_b] for elem_a in subset_a)
        # D_AB = sum(DxF[alias[elem_b], fos_idx_a] for elem_b in subset_b)
        
        ## Prepare for future execution. 
        # Compute ΣΣD(A ∪ B)
        DF[fos_idx_a_cup_b] = D_A + D_B + 2 * D_AB
        # Compute ΣD(x, A ∪ B)
        subset_a_cup_b = fos[fos_idx_a_cup_b]
        for elem_x in subset_a_cup_b
            DxF[alias[elem_x], fos_idx_a_cup_b] = DxF[alias[elem_x], fos_idx_a] + DxF[alias[elem_x], fos_idx_b]
        end

        ## Actually evaluate swaps.
        if length(subset_a) <= 2 && length(subset_b) <= 2
            # But skip the trivial cases.
            continue
        end
        for elem_a in subset_a
            for elem_b in subset_b
                # Compute new in-set 'distances'
                # Assuming D[i, i] = 0
                D_A′ = D_A + -2 * DxF[alias[elem_a], fos_idx_a] + 
                              2 * DxF[alias[elem_b], fos_idx_a] +
                             -2 * D[alias[elem_a], alias[elem_b]] 
                            # - D[alias[elem_b], alias[elem_b]] + D[alias[elem_a], alias[elem_a]]
                D_B′ = D_B + -2 * DxF[alias[elem_b], fos_idx_b] +
                              2 * DxF[alias[elem_a], fos_idx_b] +
                             -2 * D[alias[elem_a], alias[elem_b]] 
                            # - D[alias[elem_a], alias[elem_a]] + D[alias[elem_b], alias[elem_b]]
                if D_A > D_A′ && D_B > D_B′
                    # In-set distances are both smaller: swapping is an improvement.
                    # Expensive part, update values of ΣD(x, A) and ΣD(x, B) for x ∈ A ∪ B
                    for elem_x in subset_a_cup_b
                        DxF[alias[elem_x], fos_idx_a] = DxF[alias[elem_x], fos_idx_a] +
                                                        -D[alias[elem_x], alias[elem_a]] +
                                                         D[alias[elem_x], alias[elem_b]]
                        DxF[alias[elem_x], fos_idx_b] = DxF[alias[elem_x], fos_idx_b] +
                                                        -D[alias[elem_x], alias[elem_b]] +
                                                         D[alias[elem_x], alias[elem_a]]
                    end
                    # Update ΣΣD(A) and ΣΣD(B) to actually correspond to the swapped sets.
                    D_A = D_A′
                    D_B = D_B′
                    # Finally, swap the aliases of elem_a and elem_b
                    alias[elem_a], alias[elem_b] = alias[elem_b], alias[elem_a]
                    # Sanity verification, don't enable if you want to save time.
                    # @assert D_A ≈ sum(D[alias[elem_x], alias[elem_y]] for elem_x in subset_a for elem_y in subset_a)
                    # @assert D_B ≈ sum(D[alias[elem_x], alias[elem_y]] for elem_x in subset_b for elem_y in subset_b)
                    # if swap_count == 0
                    #    println("Swapping $(alias[elem_a]) with $(alias[elem_b]).\nD_A ($D_A) -> D_A′ ($D_A′)\nD_B ($D_B) -> D_B′ ($D_B′)")
                    # end
                    # Keep some statistics.
                    swap_count +=1
                end
            end
        end
    end
    # println(swap_count)
    # 'Materialize' the new fos by performing the renames in place.
    for subset in fos
        map!(f -> alias[f], subset, subset)
    end
end


##
# using BenchmarkTools

# Simply for testing!
# function create_random_distance_matrix(n :: Int64)
#     M = collect(Symmetric(rand(Float64, n, n)))
#     M[diagind(M)] .= 0
#     return M
# end

# function generate_random_test_fos(n :: Int64)
#     M = create_random_distance_matrix(n)
#     pdx=Int64[]
#     fos_set = LCP(M, UPGMA(), MersenneTwister(); parent_idx=pdx)
#     return M, [[convert(Int64, a) for a in f] for f in fos_set], findchildrenpairs(pdx))
# end


# bench = @benchmarkable rewrite_by_swap_fos(test_fos, child_pairs, M) setup=begin
#     M, test_fos, child_pairs = generate_random_test_fos(100)
# end
# run(bench)

# @benchmark generate_random_test_fos(100)
# using Distances
# using BenchmarkTools
# import Profile


# points = rand(2, 500)#[1 1 1 0 0 0; 0.2 0.3 0.4 0.2 0.3 0.4]
# D = pairwise(Euclidean(), points, dims=2)
# println(size(D, 2))
# n = size(points,2)
# es = _div64(n)+1
# FoS = sizehint!(Vector{IndexSet{es}}(undef, 0), n*2)
# ptx = Int64[]
# sizehint!(ptx, 2n)
# fos = LCP(D, UPGMA(), MersenneTwister(); parent_idx=ptx)
# tfos = [[convert(Int64, a) for a in f] for f in fos]
# child_pairs = findchildrenpairs(ptx)
# bench = @benchmarkable rewrite_by_swap_fos(copy(test_fos), child_pairs, M) setup=begin
#     M, test_fos, child_pairs = deepcopy($D), [copy(t) for t in tfos], deepcopy($child_pairs)
# end evals=1
# Profile.clear_malloc_data()
# Profile.@profile run(bench)
# ctfos = [copy(t) for t in tfos]
# rewrite_by_swap_fos(ctfos, child_pairs, D)

# bm = @profiler @benchmark LCP($(D), UPGMA(), rng; parent_idx=ptx)
# display(bm)
