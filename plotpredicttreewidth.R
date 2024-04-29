library(cowplot)
library(dplyr)
library(ggplot2)

# real data
df_real <- read.table("real_networks/attributes.csv", header=T, sep=",")
# simulated networks: file created by BDH_networks/attributes.jl
df_sim <- read.csv("BDH_networks/attributes.csv", header=T, sep=",")

# merge `df_real` and `df_sim` without removing unshared columns
df <- merge(df_real, df_sim, all=T)
df <- arrange(df, desc(Network)) # dim: 2520 x 7
# # sanity checks
# range(df$level) # 1 1835
# range(df$maxclustersize) # 3 746

set.seed(3)
# max cluster size vs no. of hybrids
p1 <- ggplot(df, aes(x = nhybrids, y = maxclustersize, color = ntips, fill = Network, shape = Network, size = Network, alpha = Network)) + 
  geom_jitter(width = 0.02, height = 0.02) +
  scale_x_log10(name = expression(paste(italic(h), ": number of hybrids"))) +
  scale_y_log10(name = "Max cluster size upper bound") +
  annotation_logticks() +
  scale_shape_manual(values = c(24, 20), guide = "none") +
  scale_size_manual(values = c(2, 1), guide = "none") +
  scale_alpha_manual(values = c(1, 0.5), guide = "none") +
  scale_fill_manual(values = c(hcl.colors(2, palette = "inferno")[2], "white"), guide = "none") +
  scale_color_binned(type = "viridis", breaks = c(2, 10, 25, 35, 50, 100, 1850),
                     direction = -1, name = expression(paste(italic(n), ": number of tips ")),
                     end = 0.8) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(color="black", fill=NA),
        legend.background=element_blank(),
        legend.key=element_rect(fill="transparent"),
        legend.box = "horizontal",
        legend.box.background=element_blank(),
        legend.position=c(0.05, 0.99),legend.justification = c("left", "top"),
        plot.title=element_text(vjust=-10, hjust=.05))
# plot(p1)

set.seed(3)
# max cluster size vs level
p2 <- ggplot(df, aes(x = level, y = maxclustersize, color = ntips, fill = Network, shape = Network, size = Network, alpha = Network)) + 
  geom_jitter(width = 0.02, height = 0.02) +
  scale_x_log10(name = "level") +
  scale_y_log10(name = "Max cluster size upper bound") +
  annotation_logticks() +
  scale_shape_manual(values = c(24, 20)) +
  scale_size_manual(values = c(2, 1)) +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_fill_manual(values = c(hcl.colors(2, palette = "inferno")[2], "white")) +
  scale_color_binned(type = "viridis", breaks = c(2, 10, 25, 35, 50, 100, 1850),
                     direction = -1, name = expression(paste(italic(n), ": number of tips ")),
                     end = 0.8, guide = "none") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(color="black", fill=NA),
        legend.background=element_blank(),
        legend.key=element_rect(color="black", fill="transparent"),
        legend.box = "horizontal",
        legend.box.background=element_blank(),
        legend.position=c(0.05, 0.99),legend.justification = c("left", "top"),
        plot.title=element_text(vjust=-10, hjust=.05))
# plot(p2)

p <- plot_grid(p1, p2, labels = c('A', 'B'), hjust = -2)
# colorblindr::cvd_grid(p) # simulate colorblindness for plot `p`

ggsave("figures/predicttreewidth.pdf", plot=p, width = 8, height = 2.5, dpi = 300)