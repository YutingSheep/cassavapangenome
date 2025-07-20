# Sequence divergence

- [Mash distance estimation (Mash)](#mash-distance-estimation-mash)
- [Mapping rate (Minimap2)](#mapping-rate-minimap2)

------

## Mash distance estimation (Mash)

Haplotype phasing was performed at the chromosome level, and sequence divergence between haplotypes was evaluated independently for each chromosome. We calculated a pairwise distance matrix using [Mash](https://mash.readthedocs.io/en/latest/) (v2.3) program. To further explore sequence divergence, we classified haplotypes into three groups: **Manihot glaziovii** (GLZ), **Manihot esculenta** ssp. **flabellifolia** (FLA), and cultivated cassava (CUL). We then extracted pairwise Mash distances from the matrix corresponding to comparisons between GLZ-CUL, FLA-CUL, and CUL-CUL to assess genetic differentiation among these groups. Mash distance analysis revealed substantial sequence divergence between **Manihot glaziovii** and cultivated cassava, which impacts node separation and edge identification in graph pangenome construction. Therefore, we excluded the two **Manihot glaziovii** (s253) haplotypes for the construction of cassava pangenome.

```shell
mashtree_bootstrap.pl --reps 100 --numcpus $threads  *.${chr}.fasta -- --min-depth 0 > mashtree.bootstrap.dnd
```

------

## Mapping rate (Minimap2)

We align 116 haplotypes (excluding the two **Manihot glaziovii** haplotypes) to the AM560 reference genome using the [Minimap2](https://github.com/lh3/minimap2) program. For each alignment file, [SAMtools](https://github.com/samtools/samtools) was used to calculate per-site coverage. The mapping rate of each haplotype was then determined by counting the detected sites divided by the reference genome size. 

```shell
# 1. Align haplotype to reference using Minimap2
# -ax asm5: assembly-to-reference mapping preset
# -t: number of threads
minimap2 -ax asm5 -t ${THREADS} ${REFERENCE} ${HAPLOTYPE} > ${OUTPUT_DIR}/alignments/${SAMPLE}.sam
    
# 2. Convert SAM to BAM and sort
# -@: number of threads
# -bS: convert SAM to BAM format
samtools view -@ ${THREADS} -bS ${OUTPUT_DIR}/alignments/${SAMPLE}.sam > ${OUTPUT_DIR}/alignments/${SAMPLE}.bam
    
# Sort BAM file by genomic coordinates
samtools sort -@ ${THREADS} -o ${OUTPUT_DIR}/alignments/${SAMPLE}.sorted.bam ${OUTPUT_DIR}/alignments/${SAMPLE}.bam
    
# 3. Calculate per-base coverage depth
# -a: output all positions including those with zero coverage
samtools depth -a ${OUTPUT_DIR}/alignments/${SAMPLE}.sorted.bam > ${OUTPUT_DIR}/coverages/${SAMPLE}.depth
    
# 4. Calculate mapping rate
# Count positions with coverage > 0 (aligned bases)
COVERED_SITES=$(awk '$3 > 0 {count++} END {print count}' ${OUTPUT_DIR}/coverages/${SAMPLE}.depth)
    
# Calculate mapping rate as: covered_sites / reference_size
# scale=4: show 4 decimal places
MAPPING_RATE=$(echo "scale=4; ${COVERED_SITES} / ${REF_SIZE}" | bc)
```

