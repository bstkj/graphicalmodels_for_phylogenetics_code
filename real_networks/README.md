# Real networks

## Data
- Networks from real data, coded in extended newick format.
The `Source` column indicates which figure of which study the network was coded
from.

|| Network | Source |
| --- | --- | --- |
|1.| `bergstrom_2020.phy` | [Fig 3a (left)](https://doi.org/10.7554/eLife.85492) |
|2.| `hajdinjak_2021.phy` | [Fig 3c (left)](https://doi.org/10.7554/eLife.85492) |
|3.| `lazaridis_2014.phy` | [Fig 3](https://doi.org/10.1038/nature13673) |
|4.| `librado_2021.phy` | [Fig 3b (left)](https://doi.org/10.7554/eLife.85492) |
|5.| `lipson_2020b.phy` | [Extended Data Fig 4](https://doi.org/10.1038/s41586-020-1929-1) |
|6.| `nielsen_2023.phy` | [Fig 3 (left)](https://doi.org/10.1371/journal.pgen.1010410) |
|7.| `muller_2022.phy` | [Fig 1a](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9297283/) |
|8.| `neureiter_2022.phy` | [Fig 5a](https://doi.org/10.1057/s41599-022-01211-7) |
|9.| `sikora_2019.phy` | [Fig 4c (left)](https://doi.org/10.7554/eLife.85492) |
|10.| `sun_2023.phy` | [Fig 4c](https://doi.org/10.1038/s41559-023-02185-8) |
|11.| `wang_2021.phy` | [Fig 4b (left)](https://doi.org/10.7554/eLife.85492) |

- `attributes.csv`: saved network attributes (e.g. upper bound for the treewidth
of the moralized network, no. of tips, no. of hybrids, level, is network
treechild)

## Scripts
- `lipson2020b_notes.jl`: documents how source network was converted to
`lipson_202b.phy` by suppressing degree-2 nodes and assigning length-0 edges a
positive edge length
- `muller2022_nexus2newick.jl`: converts source network to `muller_2022.phy` by
estimating inheritance weights from recombination breakpoints
- `attributes.jl`: creates `attributes.csv`