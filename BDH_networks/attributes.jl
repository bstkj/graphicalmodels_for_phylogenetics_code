using PhyloNetworks, PhyloGaussianBeliefProp
using DataFrames
const PGBP = PhyloGaussianBeliefProp

all_nets = readMultiTopology("BDH_networks/all_nets_nsample_10.newick")
df = DataFrame(maxclustersize=zeros(length(all_nets)),
        ntips=zeros(length(all_nets)), nhybrids=zeros(length(all_nets)))
for (i, net) in enumerate(all_nets) # Warning: this takes a long time to run
    ct = PGBP.clustergraph!(net, PGBP.Cliquetree())
    # maximum clique size ≥ treewidth(moralized network) + 1
    df[i, 1] = maximum(length(ct[lab][1]) for lab in PGBP.labels(ct))
    # extract no. of tips and no. of hybrids from network
    df[i, 2] = net.numTaxa
    df[i, 3] = net.numHybrids
end

using RCall
R"""
library(SiPhyNetwork)
library(dplyr)
df <- $df # pass `df` into R environment
# read in networks as `evonet` objects
all_nets <- lapply(1:nrow(df), function(i)
    read.net(text=scan(file="BDH_networks/all_nets_nsample_10.newick",
             what=character(), skip=i-1, nlines=1, quiet=TRUE)))
df$Network <- "simulated"
df$level <- sapply(all_nets, function(net) getNetworkLevel(net))
df$treechild <- sapply(all_nets, function(net) isTreeChild(net))
df$name <- paste0("sim_", 1:nrow(df))
# remove trees, and networks whose blobs are all simple 2-cycles
df <- df %>% filter(maxclustersize > 2) # from 3770 networks to 2509 networks

write.csv(df,file="BDH_networks/attributes.csv",row.names=F)
"""