# Translocation Identification

- [Contig community network analysis](#contig-community-network-analysis)
- [Hi-C validation of translocation events](#hi-c-validation-of-translocation-events)
- [Breakpoint analysis and homologous regions](#breakpoint-analysis-and-homologous-regions)

------

## Contig community network analysis

We performed community analysis using all contig sequences by conducting all-to-all pairwise alignments with [wfmash](https://github.com/waveygang/wfmash). To accelerate processing, we implemented parallelized paired-contig comparisons. Alignment parameters were optimized using mash-estimated haplotype distances in cassava, applying the -m option to generate approximate mappings that retained ≥95% average nucleotide identity-capturing robust signals of homologous recombination. We then calculated the total_match_ratio scores - accounting for reference/query contig lengths and aligned regions. Low-scoring alignments and contigs shorter than 1 Mb were filtered out to enhance resolution.

```shell
#Generate Alignments Using wfmash
# Use wfmash with merge (-m) to align and merge contigs  
wfmash -m reference.fa query.fa -o WFMASH.paf 
#Step 2: Convert PAF to Gephi-Compatible CSV
python calculate_paf_alignment_to_Gephi_input.py \
    -i WFMASH.paf \  
    -o Filtered_signal_for_gephi.csv  
#File Specification: Filtered_signal_for_gephi.csv
#Purpose: Input for Gephi network visualization.
#Customizable Parameters:
#Contig length filtering (adjustable via script parameters).
#Edge weights (based on alignment scores/lengths).
```

Next, we constructed an undirected graph in which nodes represented contigs, edges denoted alignments, and edge weights reflected the total_match_ratio. Chromosomal identities were assigned based on the AM560 reference genome. We visualized this network in [Gephi](https://gephi.org/) (v0.10.1)26 using the ForceAtlas2 layout algorithm.  As expected, the analysis showed that there was significant admixture of contigs (derived from different haplotypes) between the Chr01 and Chr02 communities, consistent with their extensive translocations.


------

## Hi-C validation of translocation events

s12 is a cassava accession carrying a heterozygous translocation. We mapped the Hi-C data to the haplotype which harbor the translocation.

The translocation was validated based on two key observations:

1. Intra-chromosomal Hi-C signal continuity (indicating correct haplotype assembly).
2. Inter-chromosomal Hi-C interactions at translocation breakpoints (confirming the chromosomal rearrangement).

```shell
ln -s ${HicDir}/${i}/${i}_R1.fq.gz r1Reads
ln -s ${HicDir}/${i}/${i}_R2.fq.gz r2Reads

cp /public/home/user/data_group/s12_hap1.fa contigsFasta

samtools faidx contigFasta
bwa index contigsFasta
bwa mem -5SP contigsFasta r1Reads r2Reads -t $threads | samblaster | samtools view - -@ $threads -S -h -b -F 3340 -o HiC.bam

/public/home/user/software/HapHiC/haphic pipeline contigsFasta HiC.bam 36 --threads $threads --processes $threads

/public/home/user/software/HapHiC/scripts/../utils/juicer pre -a -o out_JBAT HiC.bam ./04.build/scaffolds.agp contigsFasta.fai >out_JBAT.log 2>&1


(java -jar -Xmx240G /public/home/user/software/HapHiC/scripts/../utils/juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log | grep PRE_C_SIZE | awk '{print \$2" "\$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)
```

------

## Breakpoint analysis and homologous regions

To identify translocation breakpoints, we analyzed contigs spanning the chr01-chr02 translocation region through a multi-step refinement process. First, we aligned these contigs to the AM560 reference genome to define preliminary breakpoint intervals. For each putative translocation event, we then extracted 50 kb of flanking sequence on both sides of these initial intervals and performed targeted realignments to chr01 and chr02 using MUMmer, allowing us to narrow down the recombination breakpoint regions with higher resolution. Next, we conducted reciprocal BLASTn comparisons between the breakpoint regions of chr01 and chr02, applying a filtering strategy based on identity value. In this step, homologous mappings with higher sequence identity near the breakpoints were prioritized, as these are most likely to correspond to the homologous regions that mediated the translocation event. Finally, we identified three different highly homologous junction sequences of 2,227 bp, 2,247 bp, and 2,258 bp in length, each with more than 99% sequence identity. These homologous fragments were subsequently mapped to all translocated haplotypes, enabling the precise determination of the translocation breakpoint positions. 

```
makeblastdb -in s12HAchr01_Approximate_breakpoint_locations.fa -dbtype nucl -input_type fasta -out s12HA
blastn -query s12HAchr02_Approximate_breakpoint_locations.fa -db s12HA -out result.txt -outfmt 6 -task blastn-short -word_size 7 -evalue 0.00001 -num_threads 4
```



------

