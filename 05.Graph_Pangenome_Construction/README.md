# Graph Pangenome Construction

* [Construction of the cassava pangenome graph](#construction-of-the-cassava-pangenome-graph)

  * [Wfmash alignment](#wfmash-alignment)

  * [Seqwish graph construction](#seqwish-graph-construction)

  * [PGGB](#pggb)

* [Graph growth curve](#graph-growth-curve)

* [Non-reference nodes analysis](#non-reference-nodes-analysis)

## Construction of the cassava pangenome graph

Pangenome graph construction was performed independently for each  chromosome using the PGGB pipeline. To align the pangenome graph  construction, we manually corrected the translocation between chromosome 01 and chromosome 02 in the haplotype sequence. Specifically, all-versus-all pairwise alignments were computed using [wfmash](https://github.com/waveygang/wfmash) (v0.18.0). The resulting PAF-formatted alignments were then transformed into a graph representation in GFA format using [seqwish](https://github.com/ekg/seqwish) (v0.7.11). The partial order alignment (POA) process within [PGGB](https://github.com/pangenome/pggb) was optimized using the asm20 scoring parameters.

```shell
# Input/output setup
INPUT="haplotypes.fasta"
OUTPUT="pangenome_graph"
#one-step PGGB, including Wfmash alignment、Seqwish graph construction and PGGB
pggb -i in.pggb.fa.gz -p 90 -s 10000 -n 117 -k 47 -P 1,4,6,2,26,1 -O 0.001  -t 60 -V 'AM560' -o $PWD
#mkdir -p $OUTPUT
```

### Wfmash alignment alone

```shell
wfmash \
    -p 0.90 \          # 90% min identity
    -s 10000 \         # 10,000 bp min segment length
    -k 47 \            # k-mer size 47
    -t 8 \             # 8 threads
    $INPUT \
    > $OUTPUT/alignments.paf
```


## Graph growth curve

We characterized the growth of the cassava pangenome graph using the [gretl-GRaph Evaluation Toolkit](https://github.com/MoinSebi/gretl) (v0.1.1). To minimize sampling order effects on pangenome growth curves, we conducted a gretl bootstrap analysis with the AM560 reference genome as the initial sample, generating growth curves for both node and sequence components.

```
#Node calculation and visualization
panacus histgrowth -t4 -l 1,1,1 -q 0,0.1,0.9 -S PGGB.gfa > PGGB.node.tsv
panacus-visualize -e PGGB.node.tsv > PGGB.node.pdf
panacus histgrowth -t4 -l 1,1,1 -q 0,0.1,0.9 -S PGGB.gfa > PGGB.base.tsv
panacus-visualize -e PGGB.base.tsv > PGGB.base.pdf

#Alternative for bootstrap on each chromosome
for i in $(seq -w 1 18);do gretl bootstrap -g PGGB.chr${i}.gfa -o chr${i}.bootstrap
#Nodes/sequences with the tag ​​'S'​​ in the third column represent ​​pan-type growth values​​.
```

## Non-reference nodes analysis

For identification of non-reference sequences, we employed [odgi](https://github.com/pangenome/odgi) (v0.8.4) to detect genomic regions absent in the AM560 reference, retaining non-reference nodes larger than 50 bp. These non-reference nodes were subsequently annotated with EDTA TE classifications using the buildSummary.pl script, quantifying the proportions of major TE families including Copia/LTR, Gypsy/LTR, unknown/LTR, CACTA/TIR, Mutator/TIR, PIF_Harbinger/TIR, Tc1_Mariner/TIR, hAT/TIR, and Helitron/nonTIR elements.

```
#Extracting Path Names
odgi paths -i ${i} -L | grep 'AM560' > reference_${i}.txt
#Extracting Non-Reference Node Names​
odgi paths -i ${i} --non-reference-ranges reference_${i}.txt > non_reference_${i}.bed
​​#Correcting Chromosome Names​
sed -i -e 's/#1#/HA/g' -e 's/#2#/HB/g' non_reference_${i}.bed
#Filtering (>50 bp) to Generate BED File​
awk '{print $0, $3 - $2}' non_reference_${i}.bed > non_reference_${i}_sorted.bed
#Overlap with EDTA Annotations​​
#Example workflow:
bedtools intersect -a BRA001554_hap1.bed -b BRA001554_hap1.mod.EDTA.TEanno.gff3 -wa -wb | cut -f4-20 > BRA.tmp.gff3
perl gff2bed.pl BRA.tmp.gff3 structural > BRA.tmp.bed
perl -nle 'my ($chr, $s, $e, $anno, $dir, $supfam)=(split)[0,1,2,3,8,12]; print "10000 0.001 0.001 0.001 $chr $s $e NA $dir $anno $supfam"' BRA.tmp.bed > BRA.tmp.out
perl buildSummary.pl BRA.tmp.out > BRA.tmp.out.stat

Reference: https://github.com/oushujun/EDTA/issues/169
```

