# Graph Pangenome Construction

- [Construction of the cassava pangenome graph](#construction-of-the-cassava-pangenome-graph)
  - [Wfmash alignment](#wfmash-alignment)
  - [Seqwish graph construction](#seqwish-graph-construction)
  - [PGGB](#pggb)
- [Graph growth curve](#graph-growth-curve)
- [Non-reference nodes analysis](#non-reference-nodes-analysis)

------

## Construction of the cassava pangenome graph

Pangenome graph construction was performed independently for each  chromosome using the PGGB pipeline. To align the pangenome graph  construction, we manually corrected the translocation between chromosome 01 and chromosome 02 in the haplotype sequence. Specifically, all-versus-all pairwise alignments were computed using [wfmash](https://github.com/waveygang/wfmash) (v0.18.0). The resulting PAF-formatted alignments were then transformed into a graph representation in GFA format using [seqwish](https://github.com/ekg/seqwish) (v0.7.11). The partial order alignment (POA) process within [PGGB](https://github.com/pangenome/pggb) was optimized using the asm20 scoring parameters.

```shell
# Input/output setup
INPUT="haplotypes.fasta"
OUTPUT="pangenome_graph"
mkdir -p $OUTPUT
```

### Wfmash alignment

```shell
wfmash \
    -p 0.90 \          # 90% min identity
    -s 10000 \         # 10,000 bp min segment length
    -k 47 \            # k-mer size 47
    -t 8 \             # 8 threads
    $INPUT \
    > $OUTPUT/alignments.paf
```

### Seqwish graph construction

```shell
seqwish \
    -k 47 \            # k-mer size 47
    -t 8 \             # 8 threads
    -s $INPUT \
    -p $OUTPUT/alignments.paf \
    -g $OUTPUT/raw.gfa
```

### PGGB

```shell
pggb \
    -p $OUTPUT/final \
    -s asm20 \         # asm20 scoring preset
    -t 8 \             # 8 threads
    $OUTPUT/raw.gfa
```

------

## Graph growth curve

We characterized the growth of the cassava pangenome graph using the [gretl-GRaph Evaluation Toolkit](https://github.com/MoinSebi/gretl) (v0.1.1). To minimize sampling order effects on pangenome growth curves, we conducted a gretl bootstrap analysis with the AM560 reference genome as the initial sample, generating growth curves for both node and sequence components. 

```
#待补充
```

------

## Non-reference nodes analysis

For identification of non-reference sequences, we employed [odgi](https://github.com/pangenome/odgi) (v0.8.4) to detect genomic regions absent in the AM560 reference, retaining non-reference nodes larger than 50 bp. These non-reference nodes were subsequently annotated with EDTA TE classifications using the buildSummary.pl script, quantifying the proportions of major TE families including Copia/LTR, Gypsy/LTR, unknown/LTR, CACTA/TIR, Mutator/TIR, PIF_Harbinger/TIR, Tc1_Mariner/TIR, hAT/TIR, and Helitron/nonTIR elements.

```
#待补充
```

