#= Code to reproduce figure of network from Fig. 3 of Lazaridis et al. (2014):
Ancient human genomes suggest three ancestral populations for present-day
Europeans (https://www.nature.com/articles/nature13673)

- Major/minor hybrid edges are colored blue
- Non-root internal node labels are based on a greedy min-fill elimination order
=#
using DataFrames
using PhyloNetworks, PhyloPlots
using RCall

net = readTopology("real_networks/lazaridis_2014.phy")
# plot and save figure
R"pdf"("figures/lazaridis_2014.pdf",width=16,height=5.5)
R"par"(mar=[1.1,0,1.1,0],cex=1.2)
nodelabel = DataFrame(
    num=[-2, 19, 16, 15, 14, 13, 12, 9, 8, 7, 10, 4, 3],
    lab=["","12","11","10","8","6","5","9","7","4","2","3","1"])
plot(net,nodelabel=nodelabel,tipcex=2,nodecex=2,edgewidth=4,
    majorhybridedgecolor="blue",minorhybridedgecolor="blue",arrowlen=0.1)
R"dev.off()"

#= (1) trim out white-space bordering figure, (2) add back thinner border,
(3) save image =#
R"""
library(magick)
fig <- image_read_pdf("figures/lazaridis_2014.pdf")
image_write(image_border(image_trim(fig),"white","100x100"),
    "figures/lazaridis_2014_trim.pdf",format="pdf",density=600)
unlink("figures/lazaridis_2014.pdf") # delete untrimmed figure
"""