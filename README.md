# Leveraging graphical model techniques to study evolution on phylogenetic networks

by Benjamin Teo, Paul Bastide & Cécile Ané

This repository contains code to reproduce simulations and figures in:
"Leveraging graphical model techniques to study evolution on phylogenetic networks".

---

There are 3 top-level scripts, which plot and save (some) figures that appear in
the main text or supplementary materials:
- `plotlazaridisnet.jl`: Fig 4(a)
- `plotpredicttreewidth.R`: Fig 6
- `plotlbpaccuracy.R`: Fig 7, S2, S3

These respectively access data from `BDH_networks/`, `real_networks/` and
`lbp_accuracy/`, and save their output to `figures/`. The contents of these 4
top-level folders are summarized below, though each folder contains its own
README with more details.

- `BDH_networks/`:
    - scripts to sample simulated networks and save their attributes
    - saved network attributes
- `real_networks/`:
    - real networks and scripts to save their attributes
    - saved network attributes
    - scripts documenting how the real networks used for the "Accuracy of loopy
    BP" example were obtained from source data
- `lbp_accuracy/`:
    - scripts to run BP / loopy BP on the real networks used for the "Accuracy
    of loopy BP" example and save their estimates
    - saved estimates
- `figures/`:
    - saved figures produced by top-level scripts

## How to run the scripts

- All scripts should be run from the top-level of this repository
(i.e. `graphicalmodels_for_phylogenetics_code/`).

## Reproducing the same Julia environment

- The `Project.toml` and `Manifest.toml` files are provided so that the exact
package environment used to perform the Julia-related analyses can be
instantiated.
- To reproduce the environment:
    1. Open the Julia REPL within `graphicalmodels_for_phylogenetics_code`
    (which contains `Project.toml` and `Manifest.toml`)
    2. Enter the Pkg REPL and run `activate .`
```
(@v1.10) pkg> activate .
```