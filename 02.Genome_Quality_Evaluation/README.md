# Genome Quality Evaluation

We evaluated the quality of the assembled genomes, with a particular focus on haplotype phasing in heterozygous diploid assemblies. To systematically assess assembly quality, we categorized our evaluation into 3 Cs metrics: contiguity, completeness, and correctness. 

- [Contiguity](#contiguity)
  - [N50 & L50 (QUAST)](#n50--l50-quast)
- [Completeness](#completeness)
  - [Benchmark Universal Single Copy Orthologs (BUSCO)](#benchmark-universal-single-copy-orthologs-busco)
  - [Consensus quality (hereby noted as QV) and completeness (Merqury)](#consensus-quality-hereby-noted-as-qv-and-completeness-merqury)
- [Correctness](#correctness)
  - [Switch error and Harming error (calc_switchErr pipeline & WhatsHap)](#switch-error-and-harming-error-calcswitcherr-pipeline--whatshap)
  - [Different types of mis-assemblies (Flagger)](#different-types-of-mis-assemblies-flagger)

------

## Contiguity

### N50 & L50 (QUAST)

[QUAST](https://quast.sourceforge.net/) evaluates genome/metagenome assemblies by computing various metrics.

```shell
quast.py ${i}_haplotype.fa -t 2 -o ${i}.quast_genome
```

------

## Completeness

### Benchmark Universal Single Copy Orthologs (BUSCO)

[BUSCO](https://busco.ezlab.org/): assessing genome assembly and annotation completeness with single-copy orthologs

```shell
busco -i $GENOME/${i}_haplotype.fa -c $threads -o busco -m geno -l embryophyta_odb10 --offline
```

### Consensus quality (hereby noted as QV) and completeness (Merqury)

[Merqury](https://github.com/marbl/meryl): reference-free quality, completeness, and phasing assessment for genome assemblies.

```shell
/public/home/user/miniforge3/envs/merqury/bin/meryl k=19 count output ${i}.sr1.meryl $ILLUMINA/${i}_1_clean.fq.gz
/public/home/user/miniforge3/envs/merqury/bin/meryl k=19 count output ${i}.sr2.meryl $ILLUMINA/${i}_2_clean.fq.gz

/public/home/user/miniforge3/envs/merqury/bin/meryl union-sum output ${i}.sr.k19.meryl ${i}.sr*.meryl

/public/home/user/miniforge3/envs/merqury/bin/merqury.sh ${i}.sr.k19.meryl $GENOME/${i}_hap1.fa $GENOME/${i}_hap2.fa out_prefix
```

------

## Correctness

### Switch error and Harming error (calc_switchErr pipeline & WhatsHap)

The phasing quality in heterozygous diploids was calculated using the [calc_switchErr pipeline](https://github.com/tangerzhang/calc_switchErr#identify-signatures-bewteen-the-two-haplotypes-in-allhic-assembly) with some modifications. Specifically, we aligned HiFi reads to one of the haplotypes using [Minimap2](https://github.com/lh3/minimap2) and generated a BAM file. This BAM file was then processed with the [PEPPER-Margin-DeepVariant pipeline](https://github.com/google/deepvariant/tree/r1.8) (Mar 17, 2022) to produce a VCF file containing variant information. The resulting VCF was phased using [WhatsHap](https://whatshap.readthedocs.io/en/latest/changes.html#v2-2-2024-01-26) (v2.2) and treated as the standard VCF.  Next, we aligned another haplotype to the reference haplotype, and generated a BAM file. Variant calling was performed using bcftools. Since the identified variants were haploid, we applied a custom script to filter heterozygous sites and reformat the data to match the phased VCF (test VCF). Finally, both VCF files (one derived from HiFi reads and the other from the genome assembly) were integrated into WhatsHap program to calculate switch error and Hamming error rates.

```shell
## step1 HiFi reads构建标准数据集
mkdir 01.hifi-mapping
cd 01.hifi-mapping
cp /home/haplotype/${sample}_hap1.fa ${sample}.fasta
samtools faidx ${sample}.fasta
minimap2 -t 5 -ax map-pb --secondary=no ${sample}.fasta ${sample}.ccs.fq > ${sample}.sam
samtools view -bt ${sample}.fasta.fai ${sample}.sam > ${sample}.bam
samtools sort -@ 5 -o ${sample}.pb.sorted.bam ${sample}.bam
samtools index ${sample}.pb.sorted.bam
cd ..

## step2 两个单倍型比较得到实际组装phased 状态
mkdir 02.hap-snp-bcftools
cd 02.hap-snp-bcftools
minimap2 -a -x asm20 --cs -r2k -t 10 /home/haplotype/${sample}_hap1.fa /home/haplotype/${sample}_hap2.fa > ${sample}.hap1hap2.sam
samtools sort -O BAM -o ${sample}.hap1hap2.bam ${sample}.hap1hap2.sam
samtools index ${sample}.hap1hap2.bam
bcftools mpileup -f /home/haplotype/${sample}_hap1.fa --threads 10 -o ${sample}.hap1hap2.vcf -A -O v ${sample}.hap1hap2.bam
bcftools call -v -c ${sample}.hap1hap2.vcf -o ${sample}.hap1hap2.called.vcf
sed -i 's/1\/1/0|1/g' ${sample}.hap1hap2.called.vcf
sed -i 's/${sample}.hap1hap2.bam/${sample}/g' ${sample}.hap1hap2.called.vcf
cd ..

## step3 实际的phase和标准phase进行比较
mkdir 03.compare
cd 03.compare
#chr为染色体list
cp /home/.../chr .
cat chr| parallel -j 10 'whatshap phase --ignore-read-groups -o {}.phased.vcf --reference=../01.hifi-mapping/${sample}.fasta --chromosome ${sample}HA{} ../01.hifi-mapping/output.vcf.gz ../01.hifi-mapping/${sample}.pb.sorted.bam'
cat chr*.phased.vcf|grep -v '#'|grep PS > pb.wh.phase.vcf
cat chr01.phased.vcf | head -n 46 > ${sample}.hifi.phased.vcf
cat pb.wh.phase.vcf >> ${sample}.hifi.phased.vcf
sed -i 's/default/${sample}/g' ${sample}.hifi.phased.vcf
whatshap compare --names truth,whatshap --tsv-pairwise eval.tsv ${sample}.hifi.phased.vcf ../02.hap-snp-bcftools/${sample}.hap1hap2.called.vcf
```



### Different types of mis-assemblies (Flagger)

For a heterozygous diploid, we first merged the two haplotypes into a single genome file, constructed a 15-bp k-mer library using [Meryl](https://github.com/marbl/meryl), and generated an index file. The index file, merged genome, and HiFi reads were then input into [Winnowmap](https://github.com/marbl/Winnowmap) (v2.03) for alignment, producing an aligned BAM file. Variant detection was performed using [DeepVariant](https://github.com/google/deepvariant/tree/r1.8) (v1.5.0) to generate VCF files.  Subsequently, both the VCF and BAM files were processed using [Flagger](https://github.com/mobinasri/flagger/) (v0.4.0) to obtain a filtered BAM file. Based on this filtered BAM, we calculated the coverage of each site using SAMtools and then converted it into a coverage matrix using Flagger. The parameters were adjusted according to the sequencing depth of the HiFi reads (calculated based on each haplotype), and the coverage matrix was further transformed into a table file. The genome was divided into blocks of varying lengths and classifies them into four main categories: erroneous, duplicated, haploid (correctly assembled), and collapsed regions.

```shell
### flagger:v0.3.2
### deepvariant:v1.6.1
### Winnowmap:v2.03

### step1 align
cp /home/haplotype/${sample}_hap1.fa  ${sample}.fasta
cat /home/haplotype/${sample}_hap2.fa >> ${sample}.fasta
meryl count k=15 output merylDB ${sample}.fasta
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
samtools faidx ${sample}.fasta
winnowmap -W repetitive_k15.txt -ax map-pb -Y -L --eqx --cs ${sample}.fasta /home/HIFI/${sample}.ccs.fq.gz | samtools view -hb | samtools sort -@20 > ${sample}.sorted_1.bam

## step2 deepvariant
samtools index ${sample}.sorted_1.bam
singularity exec -B $PWD:/data deepvariant.sif /opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref /data/${sample}.fasta --reads /data/${sample}.sorted_1.bam --output_vcf /data/${sample}.flag.vcf  --make_examples_extra_args="keep_supplementary_alignments=true,min_mapping_quality=0" --num_shards=20 --dry_run=false
bcftools view -Ov -f PASS -m2 -M2 -v snps -e 'FORMAT/VAF<0.05 | FORMAT/GQ<8' ${sample}.flag.vcf > ${sample}.final.vcf

## step3 alternative alleles
filter_alt_reads -i "${sample}.sorted_1.bam" -o "${sample}_filter_alt.bam" -f "remove.bam" -v "${sample}.final.vcf" -t 20 -m 1000  -r 0.4

## step4
samtools depth -aa -Q 0 ${sample}_filter_alt.bam > read_alignment.depth
depth2cov -f ${sample}.fasta.fai -d read_alignment.depth -o read_alignment.cov
cov2counts -i read_alignment.cov -o read_alignment.counts
python3 /home/software/flagger/programs/src/fit_gmm.py --counts read_alignment.counts --cov 15  --output read_alignment.table
rm -rf result
mkdir result
find_blocks_from_table -c read_alignment.cov -t read_alignment.table  -p ./result/${sample}
awk '{sum[$1] += $3 - $2} END {for (chr in sum) print chr, sum[chr]}' ./result/${sample}.error.bed > ./result/error.stat
awk '{sum[$1] += $3 - $2} END {for (chr in sum) print chr, sum[chr]}' ./result/${sample}.collapsed.bed > ./result/collapsed.stat
awk '{sum[$1] += $3 - $2} END {for (chr in sum) print chr, sum[chr]}' ./result/${sample}.duplicated.bed > ./result/duplicated.stat
awk '{sum[$1] += $3 - $2} END {for (chr in sum) print chr, sum[chr]}' ./result/${sample}.haploid.bed > ./result/haploid.stat

python ./flagger.stat.py ./result/haploid.stat  ./result/duplicated.stat ./result/error.stat ./result/collapsed.stat --output_csv ./result/${sample}_combined_stat.csv --output_plot ./result/${sample}_chromosome_statistics_bar_chart.png
```



