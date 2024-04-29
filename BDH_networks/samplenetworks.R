# Prior to running this script, download `BDH_simulation.zip` (7.18 GB) from
# https://zenodo.org/records/8371004, unzip, and move to top-level of this
# repository (i.e. graphicalmodels_for_phylogenetics_code/).

# `BDH_simulation` has the following directory structure:
# BDH_simulation
# ├── README.txt
# ├── Supplemental_Materials.pdf
# ├── data
# └── sim_scripts

library(SiPhyNetwork)
## Find network files
all_files <- list.files(path = "BDH_simulation/data", pattern = "phys.rds", recursive = TRUE, full.names = TRUE)
length(all_files) # 377

## sample 10 networks from each scenario, write them to csv file
#  BDH_networks/all_nets_nsample_10.newick : 3770 lines
nsample <- 10
unlink(paste0("BDH_networks/all_nets_nsample_", nsample, ".newick")) # delete file if exists already
set.seed(1289)
for (ff in all_files) {
  all_nets <- readRDS(ff)
  for (net in all_nets[sample(1:length(all_nets), nsample)]) {
    write.net(net, file = paste0("BDH_networks/all_nets_nsample_", nsample, ".newick"), append = TRUE)
  }
}