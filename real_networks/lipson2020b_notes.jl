#=
Admixture graph from Extended Data Fig. 4 of
https://www.nature.com/articles/s41586-020-1929-1

We manually encoded this graph in Newick format after preprocessing it as follows:
(1) Degree-2 nodes are suppressed (e.g. a drift edge with a single child edge
that is an admixture edge is replaced by a single edge with both length (from
the drift edge) and inheritance weight (from the admixture edge))
(2) Any resulting length-0 edges were assigned length 1 (the minimum positive
edge length from the original graph)

The resulting graph was saved to the real_networks/ folder.
=#

## some checks
using PhyloNetworks
net = readTopology("real_networks/lipson_2020b.phy")
# 57 edges, 46 nodes: 12 tips, 12 hybrid nodes, 22 internal tree nodes
(e.length for e in net.edge) |> extrema # (1.0, 68.0), no length-0 edges
(e.gamma for e in net.edge if e.hybrid) |> extrema # (0.01, 0.99)
filter(n -> length(n.edge) > 3, net.node) |> length # 1 node with degree > 3
filter(n -> length(n.edge) > 3, net.node)[1].edge |> length # 4 incident edges (2 hybrid parents, 2 children)