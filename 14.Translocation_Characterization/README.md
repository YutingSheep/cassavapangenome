# Translocation Characterization

- [Estimation of translocation occurrence generation](#estimation-of-translocation-occurrence-generation)
- [Distance to centromere analysis](#distance-to-centromere-analysis)
  - [Centromeric motif presence](#centromeric-motif-presence)
  - [Repeat enrichment](#repeat-enrichment)
  - [Gene desert](#gene-desert)
- [Genetic diversity and Fst in cultivated haplotypes](#genetic-diversity-and-fst-in-cultivated-haplotypes)
- [Translocation locus linkage](#translocation-locus-linkage)

------

## Estimation of translocation occurrence generation

In the cassava translocation event, we observed a strong genetic linkage effect at the breakpoint regions, characterized by pronounced suppression of recombination. Based on analyses of newly assembled haplotype genomes, we confirmed the complete absence of recombination within these regions using two complementary approaches: (1) whole-genome alignment with Minimap2, which identified continuous and perfectly matched sequence blocks; and (2) sliding window analysis of genetic diversity patterns (detailed below). SNP variants from these recombination-suppressed regions were extracted and treated as segments with a recombination rate of r = 0, under the assumption that all polymorphisms arose from somatic mutations.
To reconstruct the evolutionary history of somatic mutations, we used the Manihot esculenta ssp. flabellifolia genotype as an ancestral reference and applied a molecular clock model to estimate mutation accumulation generation. The generation g was estimated using the following equation:

$$
g \sim \frac{N}{\mu \cdot L}
$$
where:

- N is the total number of observed mutations in the sample,
- μ is the single-nucleotide mutation rate (set to 2.7e-8 mutations/site/generation),
- L is the length of the effective genomic region analyzed. 

Translocation occurrence generation estimates were calculated separately for the chr01 and chr02 breakpoint regions. To enhance the robustness of these estimates, we applied a jackknife resampling method to construct confidence intervals. We obtained time estimates and corresponding confidence intervals per haplotype. We hypothesized that the distribution of somatic mutation ages represents the latest generation of translocation events.

```
#待补充
```



------

## Distance to centromere analysis

To evaluate the potential impact of translocation events on centromere localization, we performed a centromere prediction analysis based on the cassava AM560 reference genome. First, we scanned genome-wide potential centromeric feature sequences using the [CentIER](https://github.com/simon19891216/CentIER/tree/CentIERv2.0) (v2.0) program. Then, we conducted 21-mer frequency analysis (non-overlapping 500 kb windows, k=21) to detect highly repetitive regions. In addition, gene density was calculated within the same window size, and regions with low gene content (<5 genes per 500 kb) were flagged. By integrating these three features: centromeric motif presence, repeat enrichment, and gene deserts, we predicted centromere-associated regions of cassava chr01 and chr02. We found that the identified translocation events do not involve the core functional domains of centromeres.

### Centromeric motif presence

```shell
python CentIERv3.0.py -g Chr01.fa
python CentIERv3.0.py -g Chr02.fa
```

### Repeat enrichment

```shell
jellyfish count -m 31 -t 10 -s 1G ${split}.genome.fa
jellyfish stats mer_counts.jf
```

### Gene desert

```shell
python gene_density.py $ANNO/AM560.gff3 > AM560.500k.density
```

------

## Genetic diversity and Fst in cultivated haplotypes

We investigated genetic diversity based on a SNP dataset derived from haplotypes carrying the chr01-chr02 translocation. To ensure compliance with the VCFv4.2 format, a custom script was used to convert the detected genotype to the homozygous state. We used [popgenWindows.py](https://github.com/simonhmartin/genomics_general) to calculate genetic diversity and differentiation in non-overlapping 50-kb windows across the chr01 and chr02. As expected, we observed a reduction in genetic diversity within the chr01 inversion region and the chr01-chr02 translocation region. 

```
#待补充
```



------

## Translocation locus linkage

We conducted a GWAS based on the phased SNP profile of 54 cultivated cassava accessions, and treated the presence or absence of the translocation as a distinct binary phenotypic trait. To retain variants with high linkage disequilibrium (LD), we utilized [PLINK](https://www.cog-genomics.org/plink/) (v1.90b6.21) to construct LD blocks within 50-kb non-overlapping windows. The filtered SNPs were used for linkage analysis, employing a linear model in [GEMMA](https://github.com/genetics-statistics/GEMMA) (v0.98.5). We performed regression analysis based on this matrix with the parameter -lm 4. The resulting corrected P values were visualized using the [ggplot2](https://ggplot2.tidyverse.org/) package in R. To identify significant linkage groups, we applied the significance threshold 5×10-8.

```
#待补充
```

