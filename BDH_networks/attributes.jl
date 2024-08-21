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

# check maximum blob level for networks with h > 10
netids = findall(df.maxclustersize .> 2) # retain networks with non-trivial blobs
maxblobh_h_vec = Tuple{Int,Int}[]
for i in netids
    net = all_nets[i]
    if net.numHybrids > 10
        blobs = biconnectedComponents(net, true)
        maxblobh = 0
        for blob in blobs
            blobh = 0
            for e in blob # loop through edges in blob
                e.hybrid && (blobh += 1) # count no. of hybrid edges
            end
            # divide by 2 to get no. of hybrid nodes (assumes hybrids have 2 parents)
            (blobh > maxblobh) && (maxblobh = Int(blobh/2))
        end
        push!(maxblobh_h_vec, (maxblobh, net.numHybrids))
    end
end
diff = map(tup -> tup[2]-tup[1], maxblobh_h_vec)
length(diff) # 505 simulated networks with h > 10
sum(diff .≤ 2) # 500 out of these 505 have (no. of hybrids - max blob level) ≤ 2
maximum(diff) # (no. of hybrids - max blob level) ≤ 8