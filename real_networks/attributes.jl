using PhyloNetworks, PhyloGaussianBeliefProp
using DataFrames
const PGBP = PhyloGaussianBeliefProp

df = DataFrame(name=String[], maxclustersize=Int[], ntips=Int[], nhybrids=Int[])
for fname in filter(endswith(".phy"), readdir("real_networks/"))
    net = readTopology("real_networks/"*fname)
    ct = PGBP.clustergraph!(net, PGBP.Cliquetree()) # clique tree
    # maximum clique size ≥ treewidth(moralized network) + 1
    maxclustersize = maximum([length(ct[lab][1]) for lab in PGBP.labels(ct)])
    # extract no. of tips and no. of hybrids from network
    row = (fname, maxclustersize, net.numTaxa, net.numHybrids)
    push!(df, row)
end

using RCall
# extract network level and check if tree-child using R package `SiPhyNetwork`
R"""
library(SiPhyNetwork)
df <- $df # pass `df` into R environment
df$Network <- "empirical"
for (i in 1:nrow(df)) {
  net <- read.net(file=paste0("real_networks/", df$name[i]))
  df$level[i] <- getNetworkLevel(net)
  df$treechild[i] <- isTreeChild(net)
}

write.csv(df,file = "real_networks/attributes.csv",row.names=F)
"""