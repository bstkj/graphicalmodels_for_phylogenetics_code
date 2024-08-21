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

# check maximum blob level for networks with h > 10
# lipson_2020b.phy: h=12, muller_2022.phy: h=361, neureiter_2022.phy: h=32
for fname in ["lipson_2020b", "muller_2022", "neureiter_2022"]
  net = readTopology("real_networks/"*fname*".phy")
  blobs = biconnectedComponents(net, true)
  maxblobh = 0
  for blob in blobs
    blobh = 0
    for e in blob # loop through edges in blob
      e.hybrid && (blobh += 1) # count no. of hybrid edges
    end
    # divide by 2 to get no. of hybrid nodes (assumes hybrids have 2 parents)
    (blobh > maxblobh) && (maxblobh = blobh/2)
  end
  println("$fname, max blob level: $maxblobh, no. of hybrids: $(net.numHybrids), diff: $(net.numHybrids-maxblobh)")
end
# lipson_2020b, max blob level: 12.0, no. of hybrids: 12, diff: 0.0
# muller_2022, max blob level: 358.0, no. of hybrids: 361, diff: 3.0
# neureiter_2022, max blob level: 32.0, no. of hybrids: 32, diff: 0.0
# 2 out of 3 real networks with h > 10 have (no. of hybrids - max blob level) ≤ 2