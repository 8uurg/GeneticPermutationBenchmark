using Distributed
addprocs(4)
@everywhere include("./src/problems/BIV/biv.jl")
@everywhere include("./src/approaches/SimpleGA/IPSimpleGA.jl")
@everywhere include("./src/approaches/SimpleGA/RKSimpleGA.jl")
@everywhere insta = BIVInstance(sorted_sequential_inversion, 50)
n = 2
futures_list_4 = Vector{Future}(undef, n)
futures_list_8 = Vector{Future}(undef, n)
futures_list_16 = Vector{Future}(undef, n)
for i in 1:n
    futures_list_4[i] = @spawnat :any @timed optimize_rksimplega(bb_wrap_biv(insta), insta.n,  300.0, target_fitness=49.0, population_sizing_factor=4)
    futures_list_8[i] = @spawnat :any @timed optimize_rksimplega(bb_wrap_biv(insta), insta.n,  300.0, target_fitness=49.0, population_sizing_factor=8)
    futures_list_16[i] = @spawnat :any @timed optimize_rksimplega(bb_wrap_biv(insta), insta.n, 300.0, target_fitness=49.0, population_sizing_factor=16)
end
results_4 = fetch.(futures_list_4)
results_8 = fetch.(futures_list_8)
results_16 = fetch.(futures_list_16)

using Plots
# histogram(getindex.(results_4, 1), label="4", opacity=0.4, bar_width=1)
# histogram!(getindex.(results_8, 1), label="8", opacity=0.4, bar_width=1)
# histogram!(getindex.(results_16, 1), label="16", opacity=0.4, bar_width=1)