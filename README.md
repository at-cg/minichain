## <a name="started"></a>Getting Started

```sh
git clone https://github.com/gsc74/minichain
cd minichain && make
# Map sequence to graph
./minichain -cx lr test/MT.gfa test/MT-orangA.fa > out.gaf
```

## Table of Contents

- [Getting Started](#started)
- [Introduction](#intro)
- [Users' Guide](#uguide)
  - [Installation](#install)
  - [Sequence mapping](#map)
- [Limitations](#limit)
- [Credits](#credit)

## <a name="intro"></a>Introduction

minichain is a sequence-to-graph mapper.

## <a name="uguide"></a>Users' Guide

### <a name="install"></a>Installation

#### Dependencies
1) [gcc9][gcc9] or later version
2) [zlib][zlib]


### <a name="map"></a>Sequence mapping
minichain can be used for sequence-to-sequence mapping as well as sequence-to-graph mapping. Since minichain utilises [minigraph][minigraph] code-base, hence base level alignment is supported with [wavefront alignment algorithm][wfa].
```sh
# Map sequence to sequence
./minichain -t32 -cx lr test/MT-human.fa test/MT-orangA.fa > out.paf
# Map sequence to graph
./minichain -t32 -cx lr test/MT.gfa test/MT-orangA.fa > out.gaf
```

## <a name="limit"></a>Limitations

* Current version(v1.0) only supports acyclic [rGFA][rgfa] and [GFA][gfa1] for sequence-to-graph mapping.

* Graph Generation is not supported.

## <a name="credit"></a>Credits
minichain utilises code base of [minigraph][minigraph], which is released under MIT License.
Reference: [The design and construction of reference pangenome graphs with minigraph.][paper]

[wfa]: https://doi.org/10.1093/bioinformatics/btaa777
[paper]: https://doi.org/10.1186/s13059-020-02168-z
[minigraph]: https://github.com/lh3/minigraph
[zlib]: http://zlib.net/
[gcc9]: http://zlib.net/
[minimap2]: https://github.com/lh3/minimap2
[rgfa]: https://github.com/lh3/gfatools/blob/master/doc/rGFA.md
[gfa1]: https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md
[gaf]: https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf
[paf]: https://github.com/lh3/miniasm/blob/master/PAF.md
[gfatools]: https://github.com/lh3/gfatools
[bandage]: https://rrwick.github.io/Bandage/
[gfaviz]: https://github.com/ggonnella/gfaviz
[human-zenodo]: https://zenodo.org/record/6499594
