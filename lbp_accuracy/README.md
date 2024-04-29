# LBP accuracy

## Data
Estimates from running loopy BP on networks `lipson_2020b.phy` and
`muller_2022.phy`. The choice of cluster graph construction (e.g.
join-graph structuring vs factor graph) and regularization (e.g. R1 vs R2) are
indicated.
- `muller2022_joingraph_R1.csv`
- `muller2022_joingraph_R2.csv`
- `lipson2020b_joingraph_R1.csv`
- `lipson2020b_joingraph_R2.csv`
- `muller2022_factorgraph_R2.csv`
- `lipson2020b_factorgraph_R2.csv`

Record of each cluster's size for the join-graph structuring cluster graphs and
clique trees for networks `lipson_2020b.phy` and `muller_2022.phy`.
- `lbp_clustersizes.csv`

## Scripts
- `lbpaccuracy.jl`: creates `lbp_clustersizes.csv`,
`muller2022_joingraph_R1.csv`, `muller2022_joingraph_R2.csv`,
`lipson2020b_joingraph_R1.csv`, `lipson2020b_joingraph_R2.csv`,
`muller2022_factorgraph_R2.csv`, `lipson2020b_factorgraph_R2.csv`