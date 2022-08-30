## <a name="started"></a>Getting Started

```sh
git clone https://github.com/at-cg/minichain
cd minichain && make
# Map sequence to graph
./minichain -cx lr test/MT.gfa test/MT-orangA.fa > out.gaf
```

## Table of Contents

- [Getting Started](#started)
- [Introduction](#intro)
- [User's Guide](#uguide)
  - [Installation](#install)
  - [Sequence mapping](#map)
- [Benchmark](#benchmark)
- [Future work](#future_work)
- [Publications](#pub)

## <a name="intro"></a>Introduction

Minichain is designed to align long reads to pangenome graphs represented as DAGs. It can scale to pangenomes built from several human genome assemblies. We have incorporated a provably-good gap-sensitive co-linear chaining algorithm for filtering anchors (see [paper](#pub) for details). This algorithm enables accurate and fast long read alignments. Minichain uses seeding and base-to-base alignment code from [minigraph][minigraph].

## <a name="uguide"></a>User's Guide

### <a name="install"></a>Installation

#### Dependencies
1) [gcc9][gcc9] or later version
2) [zlib][zlib]


### <a name="map"></a>Sequence mapping
Minichain can be used for both sequence-to-sequence alignment as well as sequence-to-graph alignment. Users can run quick tests on [sample data](data/) using the following commands. The alignment output is provided in [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md) and [GAF](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf) format respectively.
```sh
# Map sequence to sequence
./minichain -cx lr test/MT-human.fa test/MT-orangA.fa > out.paf
# Map sequence to graph
./minichain -cx lr test/MT.gfa test/MT-orangA.fa > out.gaf
```

## <a name="bench"></a>Benchmark
We have compared Minichain (v1.0) with existing sequence to graph aligners to demonstrate scalability and accuracy gains. Our experiments used human pangenome graphs built by using subsets of [94 high quality haplotype assemblies](https://github.com/human-pangenomics/HPP_Year1_Assemblies) provided by the Human Pangenome Reference Consortium, and [CHM13 human genome assembly](https://www.ncbi.nlm.nih.gov/assembly/GCA_009914755.4) provided by the Telomere-to-Telomere consortium. Using a simulated long read dataset with 0.5x coverage, and graphs of three different sizes, Minichain shows superior read mapping precision. For the largest DAG from 95 haplotypes, Minichain used 24 minutes and 25 GB RAM with 32 threads.

<img src="./data/plot.png" width="700">

## <a name="future_work"></a>Future work
We plan to continue adding features in future releases. 

* Support for cyclic graphs, which can be either pangenome references or assembly graphs.

* Support for chromosome-long query sequences. This is needed for incremental graph generation.

## <a name="pub"></a>Publications

- **Ghanshyam Chandra and Chirag Jain**. "[Sequence to graph alignment using gap-sensitive co-linear chaining](https://doi.org)". *BioRxiv*, 2022.
