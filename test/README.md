# Generating haplotype-aware pangenome graphs from Human C4 Sequences

## Graph Generation
- To reproduce the results, we begin by downloading 96 Human C4 sequences from the Human Pangenome Reference Consortium (HPRC), which are hosted at [Zenodo](https://zenodo.org/records/6617246).
- We utilize the C4 sequence from the CHM13 cell line as the reference and exclude the GRCh38 cell line C4 sequence.
- For the remaining 94 sequence, we generate a pangenome graph by augmenting structural variations (SVs) of at least 50 base pairs (bp) with [Minigraph](https://github.com/lh3/minigraph).
- Subsequently, we align each sequence back to the C4 pangenome graph with [Minichain](https://github.com/at-cg/minichain) using `-cx asm` flag and augment the corresponding haplotype paths as `W` lines in the graph.
- Graphs are generated in [Graphs](Graphs) folder, [C4-CHM13_mg.gfa](C4-CHM13_mg.gfa) is rGFA pangenome graph while [C4-CHM13.gfa](C4-CHM13.gfa) is GFA v1.1 pangenome graph with `W` lines.
- Script [Gen_Graph.py](Gen_Graph.py) can be used/modified to generate pangenome graphs from Human C4 sequences stored in [Genomes](Genomes) folder.

