#=
Recombination network from Fig. 1a of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9297283/

Input: Nexus file (Source_Figure1/sars-like_all.tree) from
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9297283/bin/41467_2022_31749_MOESM6_ESM.zip
(saved to example_networks/ folder)

(1) extract network string and remove edge metadata, except for breakpoints (e.g pr={25137-29674})
(2) compute inheritance weights for recombination edges from breakpoints
        E.g. (29674-25137+1)/29675 ≈ 0.153 since breakpoints are 0-indexed, the
        corresponding major recombination edge is labeled pr={0-25136}
(3) replace breakpoints label with inheritance weights for recombination edges
        E.g. "#H4[pr={25137-29674}]:114.87" becomes "#H4:114.87::0.153"

Output: Newick file (example_networks/sars-like_all.phy)
=#

nw_str = readlines("example_networks/sars-like_all.tree")[3]
nw_str_1 = replace(nw_str,r"tree STATE_0 = (?<nw>.*)" => s"\g<nw>")
nw_str_2 = replace(nw_str_1,r"&loci=\{[^}]*\}" => "") # remove &loci={..., ...}
nw_str_3 = replace(nw_str_2,r"length=[\d.E-]*" => "") # remove length=...
nw_str_4 = replace(nw_str_3,r"posterior=[\d.E-]*" => "") # remove posterior=...
nw_str_5 = replace(nw_str_4,r"height_95%_HPD=\{[^}]*\}" => "") # remove height_95%_HPD={..., ...}
nw_str_6 = replace(nw_str_5,r"\[[,]+\]" => "") # remove any commas in [...]
nw_str_7 = replace(nw_str_6,r"\[[,pr={]*(?<range>[\d-]+)[,}]*\]" => s"[\g<range>]") # extract range from pr={...}

subst = [] # substitutions to annotate network string with inheritance weights
for (i, rmatch) in enumerate(eachmatch(r"#H(?<hybridnum>[\d]+)\[(?<breakpoints>[\d-]+)\]:(?<edgelength>[\d.]+)", nw_str_7))
    hybridnum = rmatch.captures[1]
    breakpoints = rmatch.captures[2]
    edgelength = rmatch.captures[3]
    gamma = round(-(eval(Meta.parse(breakpoints))-1)/29675,digits=3)
    push!(subst, "#H$hybridnum[$breakpoints]:$edgelength" => "#H$hybridnum:$edgelength::$gamma")
end
nw_str_8 = replace(nw_str_7, subst...)

file_path = "example_networks/muller_2022.phy"
open(file_path, "w") do file
    write(file, nw_str_8)
end
# manually removed the length-0 root edge (e.g. "(...):0.0;" to "(...);"

## some checks
using PhyloNetworks
net = readTopology("example_networks/muller_2022.phy")
# 1161 edges, 801 nodes: 40 tips, 361 hybrid nodes, 400 internal tree nodes
(e.length for e in net.edge) |> extrema # (0.003169596180669032, 907.3366274607556), no length-0 edges
filter(e -> e.hybrid && e.isMajor && e.gamma < 0.5, net.edge) # [], major hybrid edges assigned by gamma
filter(n -> length(n.edge) > 3, net.node) # [], no polytomies