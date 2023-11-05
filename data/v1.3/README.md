# Reproducing Results for Minichain v1.3

## Overview of the Scripts for Result Reproduction

**Step 1:** Data Preparation
- To reproduce the results, we begin by downloading 61 major histocompatibility complex (MHC) genomes from the Human Pangenome Reference Consortium (HPRC), which are available on [Zenodo](https://zenodo.org/records/6617246).
- We utilize the MHC genome from the CHM13 cell line as the reference and exclude the GRCh38 cell line MHC genome.
- For the remaining 59 haplotypes, we generate a pangenome graph by augmenting structural variations (SVs) of at least 50 base pairs (bp) with [Minigraph](https://github.com/lh3/minigraph).
- Subsequently, we align each haplotype back to the pangenome graph with [Minigraph](https://github.com/lh3/minigraph) using `-cx asm` flag and augment the corresponding haplotype paths as `W` lines in the graph.

**Step 2:** Simulating Haplotypes
- We simulate 135 MHC haplotypes, creating a mosaic of haplotype paths derived from the graph and recording the true haplotype paths.
- We split the simulated haplotypes into three groups, each consisting of 45 haplotypes, with substitution error rates of 0.1%, 1%, and 5%.

**Step 3:** Alignment and Haplotype Switches
- The simulated MHC haplotypes are aligned to the MHC pangenome graph, and the chained haplotype paths are recorded.

**Step 4:** Evaluation Metrics
- Using the information from true haplotype paths and chained haplotype paths, we calculate both the Pearson correlation coefficient and F1-scores.

**Step 5:** Result Visualization
- The generated plots are saved in a folder labeled "align_x," where 'x' represents the substitution rate.

To reproduce the results, execute the following commands:

```bash
chmod +x Reproduce.py
./Reproduce.py -t4
```

Here, `-t4` represents the thread count, and memory requirements increase proportionally at a rate of `15` GB per thread.

