# Reproducing results for Minichain v1.3
We assume that [conda](https://docs.conda.io/en/latest/) is installed and available in your path. To confirm, please run `which conda` 
## To reproduce the results, execute the following commands:
```bash
# Inside minichain/data/v1.3 folder
chmod +x Reproduce.py
./Reproduce.py -t8
```
Here `-t8` represents the thread count. The memory required is about `15` GB per thread. We recommend  setting the maximum thread count based on your hardware.

## Benchmarking methodolody

**Step 1:** Data preparation
- Download 61 major histocompatibility complex (MHC) sequences from the Human Pangenome Reference Consortium (HPRC). These are available on [Zenodo](https://zenodo.org/records/6617246).
- We use the MHC sequence from CHM13 T2T assembly as the reference during graph construction. We do not use the MHC sequence from GRCh38.
- Using 60 haplotypes, we build a pangenome graph by augmenting structural variations (SVs) with [Minigraph](https://github.com/lh3/minigraph).
- Subsequently, we align each haplotype back to the pangenome graph with [Minigraph](https://github.com/lh3/minigraph) using `-cx asm` flag. We use these alignments to store haplotype paths in the GFA file.

**Step 2:** Simulating haplotypes
- We simulate 135 MHC haplotype sequences. Each sequence is a mosaic of the haplotype paths stored in the graph.
- We add substitution errors with rates 0.1%, 1%, and 5% in the query sequences.

**Step 3:** Alignment and haplotype switches
- The simulated MHC haplotypes are aligned to the MHC pangenome graph, and the chained haplotype paths are recorded with `-b1` parameter.

**Step 4:** Evaluation metrics
- Using the known ground truth (recombination events) and the chained haplotype paths, we measure Pearson correlation coefficient and F1-scores.

**Step 5:** Result visualization
- The generated plots will be saved in a folder labeled "align_x," where 'x' represents the substitution rate.

