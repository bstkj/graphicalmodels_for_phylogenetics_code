# BDH networks

## Data
- `all_nets_nsample_10.newick`: sampled simulated birth-death-hybridization
(BDH) networks (10 per simulation setting) from [*Exploring the Distribution of
Phylogenetic Networks Generated Under a Birth-Death-Hybridization
Process*](https://doi.org/10.18061/bssb.v2i3.9285)
- `attributes.csv`: saved network attributes (e.g. upper bound for the
treewidth of the moralized network, no. of tips, no. of hybrids, level, is
network treechild)

## Scripts
- `samplenetworks.R`: creates `all_nets_nsample_10.newick`
- `attributes.jl`: creates `attributes.csv`