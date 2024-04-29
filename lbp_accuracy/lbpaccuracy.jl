# Reproduce data for "Accuracy of loopy BP" Fig 7, S2, S3

using BenchmarkTools
using CSV, DataFrames
using PhyloGaussianBeliefProp, PhyloNetworks
using MetaGraphsNext
using LinearAlgebra
using StatsBase
using Random
using Tables
const PGBP = PhyloGaussianBeliefProp

################################################################################
#= Network from Fig 1a of Müller et al. (2022): A Bayesian approach to infer
recombination patterns in coronaviruses =#

## simulate tip data
net = readTopology("real_networks/muller_2022.phy");
preorder!(net);
Random.seed!(3);
sim = simulate(net, ParamsBM(0, 1)); # BM model with parameters μ=0, σ²=1
df = DataFrame(taxon=tipLabels(net), x=sim[:Tips]);
tbl_x = columntable(select(df, :x));

## clique tree
ct = PGBP.clustergraph!(net, PGBP.Cliquetree());
summarystats(map(cl -> length(ct[cl][1]), labels(ct))) # summary of cluster sizes
# Mean:6.73, Std. Dev:6.12, Min:2, Q1:4, Median:5, Q3:7, Max:54
m = PGBP.UnivariateBrownianMotion(1, 0, Inf); # 𝒩(0,∞) prior on root mean μ
b = PGBP.init_beliefs_allocate(tbl_x, df.taxon, net, ct, m);
PGBP.init_beliefs_assignfactors!(b, m, tbl_x, df.taxon, net.nodes_changed);
ctb = PGBP.ClusterGraphBelief(b);
spt = PGBP.spanningtree_clusterlist(ct, net.nodes_changed);
PGBP.calibrate!(ctb, [spt]);
# Benchmarked on 02/14/24:
# BenchmarkTools.DEFAULT_PARAMETERS.samples = 20
# BenchmarkTools.DEFAULT_PARAMETERS.seconds = 100
# @benchmark PGBP.calibrate!($ctb, [$spt], 100; auto=false, info=false)
# BenchmarkTools.Trial: 20 samples with 1 evaluation.
#  Range (min … max):  1.336 s …   1.562 s  ┊ GC (min … max): 0.00% … 9.49%
#  Time  (median):     1.470 s              ┊ GC (median):    7.58%
#  Time  (mean ± σ):   1.460 s ± 59.526 ms  ┊ GC (mean ± σ):  6.82% ± 3.72%

#   ▁   ▁ ▁           ▁     ▁   ▁ ▁▁▁▁█ ▁▁▁ ▁█              █  
#   █▁▁▁█▁█▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁█▁▁▁█▁█████▁███▁██▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
#   1.34 s         Histogram: frequency by time        1.56 s <

#  Memory estimate: 971.53 MiB, allocs estimate: 13080601.
# 1.460 / (length(spt[1]) * 2 * 100) ≈ 11.0 μs per message
root_ind = findfirst(b -> 1 ∈ PGBP.nodelabels(b), b); # 258
(b[root_ind].J \ I)[end,end] # root conditional variance: 1425.651004065571
PGBP.integratebelief!(b[root_ind])[1][end] # root conditional mean: -15.696383011725299
PGBP.integratebelief!(b[root_ind])[2] # log normalization constant (root cluster) == llscore: -85.04744054045547
PGBP.factored_energy(ctb)[3] # -85.0474405404234

## join-graph structuring
cg = PGBP.clustergraph!(net, PGBP.JoinGraphStructuring(10));
summarystats(map(cl -> length(cg[cl][1]), labels(cg))) # summary of cluster sizes
# Mean:6.04, Std. Dev:2.18, Min:1, Q1:4, Median:6, Q3:8, Max:10
m = PGBP.UnivariateBrownianMotion(1, 0, Inf); # 𝒩(0,∞) prior on root mean μ
b_jg = PGBP.init_beliefs_allocate(tbl_x, df.taxon, net, cg, m);
PGBP.init_beliefs_assignfactors!(b_jg, m, tbl_x, df.taxon, net.nodes_changed);
cgb = PGBP.ClusterGraphBelief(b_jg);
sch = PGBP.spanningtrees_clusterlist(cg, net.nodes_changed);
# PGBP.regularizebeliefs_onschedule!(cgb, cg);
# Benchmarked on 02/14/24:
# BenchmarkTools.DEFAULT_PARAMETERS.samples = 20
# BenchmarkTools.DEFAULT_PARAMETERS.seconds = 100
# @benchmark PGBP.calibrate!($cgb, $sch, 100; auto=false, info=false)
# BenchmarkTools.Trial: 20 samples with 1 evaluation.
#  Range (min … max):  4.829 s …   5.174 s  ┊ GC (min … max): 3.19% … 6.02%
#  Time  (median):     4.978 s              ┊ GC (median):    4.66%
#  Time  (mean ± σ):   4.989 s ± 90.683 ms  ┊ GC (mean ± σ):  4.71% ± 1.00%

#   ▁      ▁  ▁    ▁▁█  ▁ ▁▁  ▁▁  █▁▁         ▁      ▁  ▁   ▁  
#   █▁▁▁▁▁▁█▁▁█▁▁▁▁███▁▁█▁██▁▁██▁▁███▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁█▁▁█▁▁▁█ ▁
#   4.83 s         Histogram: frequency by time        5.17 s <

#  Memory estimate: 3.18 GiB, allocs estimate: 53344200.
# 4.989 / (length(sch)*length(sch[1][1]) * 2 * 100) ≈ 8.32 μs per message

## save data for plotting
### join-graph structuring + regularization algorithm R1
PGBP.init_beliefs_assignfactors!(b_jg, m, tbl_x, df.taxon, net.nodes_changed);
PGBP.regularizebeliefs_bynodesubtree!(cgb, cg);
root_ind = findfirst(b -> 1 ∈ PGBP.nodelabels(b), b_jg); # 98
res = DataFrame(posvar=Float64[], posmean=Float64[], norm=Float64[], fenergy=Float64[]);
for i in 1:400
    PGBP.calibrate!(cgb, sch)
    row = (try (b_jg[root_ind].J \ I)[end,end] catch ex NaN end,
        try PGBP.integratebelief!(b_jg[root_ind])[1][end] catch ex NaN end,
        try PGBP.integratebelief!(b_jg[root_ind])[2] catch ex NaN end,
        try PGBP.factored_energy(cgb)[3] catch ex NaN end)
    push!(res, row)
end
CSV.write("lbp_accuracy/muller2022_joingraph_R1.csv", res)

### join-graph structuring + regularization algorithm R2
PGBP.init_beliefs_assignfactors!(b_jg, m, tbl_x, df.taxon, net.nodes_changed);
PGBP.regularizebeliefs_onschedule!(cgb, cg);
root_ind = findfirst(b -> 1 ∈ PGBP.nodelabels(b), b_jg); # 98
res = DataFrame(posvar=Float64[], posmean=Float64[], norm=Float64[], fenergy=Float64[]);
for i in 1:400
    PGBP.calibrate!(cgb, sch)
    row = ((b_jg[root_ind].J \ I)[end,end],
        PGBP.integratebelief!(b_jg[root_ind])[1][end],
        PGBP.integratebelief!(b_jg[root_ind])[2],
        PGBP.factored_energy(cgb)[3])
    push!(res, row)
end
CSV.write("lbp_accuracy/muller2022_joingraph_R2.csv", res)

## bethe cluster graph / factor graph
cg = PGBP.clustergraph!(net, PGBP.Bethe());
summarystats(map(cl -> length(cg[cl][1]), labels(cg))) # summary of cluster sizes
# Mean:1.74, Std. Dev:0.81, Min:1, Q1:1, Median:2, Q3:2, Max:3
m = PGBP.UnivariateBrownianMotion(1, 0, Inf); # 𝒩(0,∞) prior on root mean μ
b_fg = PGBP.init_beliefs_allocate(tbl_x, df.taxon, net, cg, m);
PGBP.init_beliefs_assignfactors!(b_fg, m, tbl_x, df.taxon, net.nodes_changed);
cgb = PGBP.ClusterGraphBelief(b_fg);
sch = PGBP.spanningtrees_clusterlist(cg, net.nodes_changed);
# PGBP.regularizebeliefs_onschedule!(cgb, cg);
# Benchmarked on 02/28/24:
# BenchmarkTools.DEFAULT_PARAMETERS.samples = 20
# BenchmarkTools.DEFAULT_PARAMETERS.seconds = 100
# @benchmark PGBP.calibrate!($cgb, $sch, 100; auto=false, info=false)
# BenchmarkTools.Trial: 20 samples with 1 evaluation.
#  Range (min … max):  2.757 s …   2.806 s  ┊ GC (min … max): 6.89% … 6.92%
#  Time  (median):     2.768 s              ┊ GC (median):    6.88%
#  Time  (mean ± σ):   2.775 s ± 15.120 ms  ┊ GC (mean ± σ):  6.87% ± 0.04%

#         █                                                    
#   ▇▁▇▁▁▁█▇▇▇▇▇▇▁▁▁▁▁▁▁▁▁▇▇▁▁▁▁▁▁▇▁▇▁▁▇▁▁▁▁▁▇▇▁▁▁▁▁▁▁▁▁▁▇▁▁▇ ▁
#   2.76 s         Histogram: frequency by time        2.81 s <

#  Memory estimate: 2.31 GiB, allocs estimate: 42076000.
# 2.775 / (length(sch)*length(sch[1][1]) * 2 * 100) ≈ 4.46 μs per message

## save data for plotting
### factor graph + regularization algorithm R2
PGBP.init_beliefs_assignfactors!(b_fg, m, tbl_x, df.taxon, net.nodes_changed);
PGBP.regularizebeliefs_onschedule!(cgb, cg);
root_ind = findfirst(b -> 1 ∈ PGBP.nodelabels(b), b_fg); # 797
res = DataFrame(posvar=Float64[], posmean=Float64[], norm=Float64[], fenergy=Float64[]);
for i in 1:400
    PGBP.calibrate!(cgb, sch)
    row = (try (b_fg[root_ind].J \ I)[end,end] catch ex NaN end,
        try PGBP.integratebelief!(b_fg[root_ind])[1][end] catch ex NaN end,
        try PGBP.integratebelief!(b_fg[root_ind])[2] catch ex NaN end,
        try PGBP.factored_energy(cgb)[3] catch ex NaN end)
    push!(res, row)
end
CSV.write("lbp_accuracy/muller2022_factorgraph_R2.csv", res)

################################################################################
#= Network from Extended Fig 4 of Lipson et al. (2020b): Ancient West African
foragers in the context of African population history =#

## simulate tip data
net = readTopology("real_networks/lipson_2020b.phy");
preorder!(net);
Random.seed!(3);
sim = simulate(net, ParamsBM(0, 1)); # BM model with parameters μ=0, σ²=1
df = DataFrame(taxon=tipLabels(net), x=sim[:Tips]);
tbl_x = columntable(select(df, :x));

## clique tree
ct = PGBP.clustergraph!(net, PGBP.Cliquetree());
summarystats(map(cl -> length(ct[cl][1]), labels(ct))) # summary of cluster sizes
# Mean:3.45, Std. Dev:1.45, Min:2, Q1:2, Median:3, Q3:4, Max:7
m = PGBP.UnivariateBrownianMotion(1, 0, Inf); # 𝒩(0,∞) prior on root mean μ
b = PGBP.init_beliefs_allocate(tbl_x, df.taxon, net, ct, m);
PGBP.init_beliefs_assignfactors!(b, m, tbl_x, df.taxon, net.nodes_changed);
ctb = PGBP.ClusterGraphBelief(b);
spt = PGBP.spanningtree_clusterlist(ct, net.nodes_changed);
PGBP.calibrate!(ctb, [spt]);
# Benchmarked on 02/14/24:
# BenchmarkTools.DEFAULT_PARAMETERS.samples = 20
# BenchmarkTools.DEFAULT_PARAMETERS.seconds = 100
# @benchmark PGBP.calibrate!($ctb, [$spt], 100; auto=false, info=false)
# BenchmarkTools.Trial: 20 samples with 1 evaluation.
# Range (min … max):  45.140 ms …  47.109 ms  ┊ GC (min … max): 0.00% … 0.00%
# Time  (median):     45.270 ms               ┊ GC (median):    0.00%
# Time  (mean ± σ):   45.412 ms ± 437.725 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

#  ▁█▁▁▁     ▁                                                   
#  █████▆▆▁▁▆█▁▁▁▆▁▁▁▁▁▁▁▆▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▆ ▁
#  45.1 ms         Histogram: frequency by time         47.1 ms <

# Memory estimate: 34.13 MiB, allocs estimate: 578801.
# 45.412 / 1000 / (length(spt[1]) * 2 * 100) ≈ 5.82 μs per message
root_ind = findfirst(b -> 1 ∈ PGBP.nodelabels(b), b); # 33
(b[root_ind].J \ I)[end,end] # root conditional variance: 31.032000074578907
PGBP.integratebelief!(b[root_ind])[1][end] # root conditional mean: 1.7010205811114112
PGBP.integratebelief!(b[root_ind])[2] # log normalization constant (root cluster) == llscore: -35.89254355788221
PGBP.factored_energy(ctb)[3] # -35.892543557882675

## join-graph structuring
cg = PGBP.clustergraph!(net, PGBP.JoinGraphStructuring(3));
summarystats(map(cl -> length(cg[cl][1]), labels(cg)))
# Mean:2.43, Std. Dev:0.53, Min:1, Q1:2, Median:2, Q3:3, Max:3
m = PGBP.UnivariateBrownianMotion(1, 0, Inf); # 𝒩(0,∞) prior on root mean μ
b_jg = PGBP.init_beliefs_allocate(tbl_x, df.taxon, net, cg, m);
PGBP.init_beliefs_assignfactors!(b_jg, m, tbl_x, df.taxon, net.nodes_changed);
cgb = PGBP.ClusterGraphBelief(b_jg);
sch = PGBP.spanningtrees_clusterlist(cg, net.nodes_changed);
# PGBP.regularizebeliefs_onschedule!(cgb, cg);
# Benchmarked on 02/14/24:
# BenchmarkTools.DEFAULT_PARAMETERS.samples = 20
# BenchmarkTools.DEFAULT_PARAMETERS.seconds = 100
# @benchmark PGBP.calibrate!($cgb, $sch, 100; auto=false, info=false)
# BenchmarkTools.Trial: 20 samples with 1 evaluation.
#  Range (min … max):  101.212 ms … 181.806 ms  ┊ GC (min … max): 0.00% … 42.31%
#  Time  (median):     102.746 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   107.033 ms ±  17.652 ms  ┊ GC (mean ± σ):  3.59% ±  9.46%

#    █                                                             
#   ▇█▆▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃ ▁
#   101 ms           Histogram: frequency by time          182 ms <

#  Memory estimate: 88.77 MiB, allocs estimate: 1558400.
# 107.033 / 1000 / (length(sch)*length(sch[1][1]) * 2 * 100) ≈ 4.87 μs per message

## save data for plotting
### join-graph structuring + regularization algorithm R1
PGBP.init_beliefs_assignfactors!(b_jg, m, tbl_x, df.taxon, net.nodes_changed);
PGBP.regularizebeliefs_bynodesubtree!(cgb, cg);
root_ind = findfirst(b -> 1 ∈ PGBP.nodelabels(b), b_jg); # 25
res = DataFrame(posvar=Float64[], posmean=Float64[], norm=Float64[], fenergy=Float64[]);
for i in 1:20
    PGBP.calibrate!(cgb, sch)
    row = ((b_jg[root_ind].J \ I)[end,end],
        PGBP.integratebelief!(b_jg[root_ind])[1][end],
        PGBP.integratebelief!(b_jg[root_ind])[2],
        PGBP.factored_energy(cgb)[3])
    push!(res, row)
end
CSV.write("lbp_accuracy/lipson2020b_joingraph_R1.csv", res)

### join-graph structuring + regularization algorithm R2
PGBP.init_beliefs_assignfactors!(b_jg, m, tbl_x, df.taxon, net.nodes_changed);
PGBP.regularizebeliefs_onschedule!(cgb, cg);
root_ind = findfirst(b -> 1 ∈ PGBP.nodelabels(b), b_jg); # 25
res = DataFrame(posvar=Float64[], posmean=Float64[], norm=Float64[], fenergy=Float64[]);
for i in 1:20
    PGBP.calibrate!(cgb, sch)
    row = ((b_jg[root_ind].J \ I)[end,end],
        PGBP.integratebelief!(b_jg[root_ind])[1][end],
        PGBP.integratebelief!(b_jg[root_ind])[2],
        PGBP.factored_energy(cgb)[3])
    push!(res, row)
end
CSV.write("lbp_accuracy/lipson2020b_joingraph_R2.csv", res)

## bethe cluster graph / factor graph
cg = PGBP.clustergraph!(net, PGBP.Bethe());
summarystats(map(cl -> length(cg[cl][1]), labels(cg))) # summary of cluster sizes
# Mean:1.72, Std. Dev:0.71, Min:1, Q1:1, Median:2, Q3:2, Max:3
m = PGBP.UnivariateBrownianMotion(1, 0, Inf); # 𝒩(0,∞) prior on root mean μ
b_fg = PGBP.init_beliefs_allocate(tbl_x, df.taxon, net, cg, m);
PGBP.init_beliefs_assignfactors!(b_fg, m, tbl_x, df.taxon, net.nodes_changed);
cgb = PGBP.ClusterGraphBelief(b_fg);
sch = PGBP.spanningtrees_clusterlist(cg, net.nodes_changed);
# PGBP.regularizebeliefs_onschedule!(cgb, cg);
# Benchmarked on 02/28/24:
# BenchmarkTools.DEFAULT_PARAMETERS.samples = 20
# BenchmarkTools.DEFAULT_PARAMETERS.seconds = 100
# @benchmark PGBP.calibrate!($cgb, $sch, 100; auto=false, info=false)
# BenchmarkTools.Trial: 20 samples with 1 evaluation.
#  Range (min … max):  133.886 ms … 204.104 ms  ┊ GC (min … max): 5.81% … 8.33%
#  Time  (median):     141.092 ms               ┊ GC (median):    8.13%
#  Time  (mean ± σ):   145.150 ms ±  15.607 ms  ┊ GC (mean ± σ):  7.17% ± 1.45%

#     ▂   █                                                        
#   ▅▅██▅▅█▁▁▅▁█▁▁▅▁▁▁▁▁▁▁▁▁▁▁▁▁▅▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▅ ▁
#   134 ms           Histogram: frequency by time          204 ms <

#  Memory estimate: 114.85 MiB, allocs estimate: 2030400.
# 145.150 / 1000 / (length(sch)*length(sch[1][1]) * 2 * 100) ≈ 4.65 μs per message

## save data for plotting
### factor graph + regularization algorithm R2
PGBP.init_beliefs_assignfactors!(b_fg, m, tbl_x, df.taxon, net.nodes_changed);
PGBP.regularizebeliefs_onschedule!(cgb, cg);
root_ind = findfirst(b -> 1 ∈ PGBP.nodelabels(b), b_fg); # 44
res = DataFrame(posvar=Float64[], posmean=Float64[], norm=Float64[], fenergy=Float64[]);
for i in 1:20
    PGBP.calibrate!(cgb, sch)
    row = ((b_fg[root_ind].J \ I)[end,end],
        PGBP.integratebelief!(b_fg[root_ind])[1][end],
        PGBP.integratebelief!(b_fg[root_ind])[2],
        PGBP.factored_energy(cgb)[3])
    push!(res, row)
end
CSV.write("lbp_accuracy/lipson2020b_factorgraph_R2.csv", res)

################################################################################
## Save cluster sizes (for boxplots)
net = readTopology("real_networks/muller_2022.phy");
ct = PGBP.clustergraph!(net, PGBP.Cliquetree());
muller_ct = map(cl -> length(ct[cl][1]), labels(ct));
cg = PGBP.clustergraph!(net, PGBP.JoinGraphStructuring(10));
muller_cg = map(cl -> length(cg[cl][1]), labels(cg));

net = readTopology("real_networks/lipson_2020b.phy");
ct = PGBP.clustergraph!(net, PGBP.Cliquetree());
lipson_ct = map(cl -> length(ct[cl][1]), labels(ct));
cg = PGBP.clustergraph!(net, PGBP.JoinGraphStructuring(3));
lipson_cg = map(cl -> length(cg[cl][1]), labels(cg));

CSV.write("lbp_accuracy/lbpaccuracy_clustersizes.csv",
    DataFrame(clustersize=vcat(muller_ct, muller_cg, lipson_ct, lipson_cg),
        grp=vcat(repeat([1],length(muller_ct)), repeat([2],length(muller_cg)),
            repeat([3],length(lipson_ct)), repeat([4],length(lipson_cg)))))