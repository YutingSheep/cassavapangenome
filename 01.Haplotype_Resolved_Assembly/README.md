# Haplotype-Resolved Assembly

Haplotype-resolved assemblies for 57 cassava accessions were generated using HiFi and Hi-C reads with [Hifiasm](https://github.com/chhylp123/hifiasm) (v0.19.8) under default parameters. Contigs derived from plasmid sequences or composed mainly of microsatellite repeats were filtered. Scaffolding was performed independently for each haplotype through two rounds of iterations to ensure the completeness and accuracy of the assembly. Hi-C reads were aligned to each haplotype using [Chromap](https://github.com/haowenz/chromap) (v0.2.5) with default settings, and contigs were ordered and oriented into 18 chromosome-scale pseudomolecules using [YaHS](https://github.com/c-zhou/yahs) (v1.1). Hi-C contact maps were visualized and manually curated to correct assembly errors and refine contig order using [Juicebox](https://github.com/aidenlab/Juicebox) (v2.20.00). To further validate phasing and scaffolding, Hi-C reads were aligned to the merged haplotypes, retaining high-quality contact signals using [HapHiC](https://github.com/zengxiaofei/HapHiC) (v1.0.2). For the reference genome of the third-generation inbred line AM560, primary contigs were assembled from ONT ultra-long reads, HiFi reads, and Hi-C data using Hifiasm and default parameters, followed by Hi-C-based scaffolding as described above.

* [Assembly](#assembly)

  * [Reference genome assembly](#reference-genome-assembly)

  * [Haplotype-resolved assembly](#haplotype-resolved-assembly)

* [Contig sequence filtering](#contig-sequence-filtering)

* [Haplotype-resolved scaffolding](#haplotype-resolved-scaffolding)

  * [Scaffolding two phased contigs separately](#scaffolding-two-phased-contigs-separately)

  * [Combine the separately scaffolded haplotypes and perform Hi-C data realignment](#combine-the-separately-scaffolded-haplotypes-and-perform-hi-c-data-realignment)

  * [Manual scaffold adjustment](#manual-scaffold-adjustment)

* [Assembly Result Visualization](#assembly-result-visualization)

## Assembly

### Reference genome assembly

```shell
hifiasm -t $threads -o TEST1-MS --primary -l3 --h1 ${Hic_1} --h2 ${Hic_2} --ul $ONT $CCS
```

### Haplotype-resolved assembly

```shell
hifiasm -t $threads -o TEST1-MS --primary -l3 --h1 ${Hic_1} --h2 ${Hic_2}  $CCS
```

## Contig sequence filtering

Since the initial contig assembly contained plastid genome sequences, we downloaded previously published cassava mitochondrial (MK176513.1) and chloroplast (EU117376.1) genome sequences from NCBI (GCF_001659605.2) and annotated our assembled contigs based on sequence alignment using Minimap2. Contigs predominantly composed of plastid DNA were excluded. Satellite repeats were mainly identified using [Meryl](https://github.com/marbl/meryl) (v1.4.1) with 41-bp K-mers. GC content analysis was also calculated in identifying satellite contigs, as well as contigs with potential assembly algorithm bias or sequencing bias. After these filtering steps, the remaining contigs were used for downstream Hi-C scaffolding.

```shell
#need minimap2 and meryl
#The parameters are optimized for cassava genome
python contig_filter.py -f draft.contig.fa -c chloroplast.fa -m mitochondrial.fa -t 10 -o contig.fa
```

## Haplotype-resolved scaffolding

### Scaffolding two phased contigs separately

```shell
ln -s /home/${sample}.clean.hap1.CLEAN.fa contigsFasta
#ln -s /home/${sample}.clean.hap2.CLEAN.fa contigsFasta

ln -s /home/HIC/${sample}_R1.fq.gz r1Reads
ln -s /home/HIC/${sample}_R2.fq.gz r2Reads

#indec
samtools faidx contigsFasta
chromap -i -r contigsFasta -o contigs.index

# alignment
chromap --preset hic -r contigsFasta -x contigs.index --remove-pcr-duplicates -1 r1Reads -2 r2Reads --SAM -o aligned.sam -t 20

#sort  
samtools view -bh aligned.sam | chromapyahs/bin/samtools sort -@ 20 -n > aligned.bam
rm aligned.sam

#step2: scaffolding
yahs contigsFasta aligned.bam
juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final.agp contigsFasta.fai
asm_size=$(awk '{s+=$2} END{print s}' contigsFasta.fai)
java -Xmx36G -jar juicer_tools_1.19.02.jar pre out_JBAT.txt out_JBAT.hic <(echo "assembly ${asm_size}")
```

### Combine the separately scaffolded haplotypes and perform Hi-C data realignment

```shell
ln -s /home/HIC/${sample}_R1.fq.gz r1Reads
ln -s /home/HIC/${sample}_R2.fq.gz r2Reads
cp /home/haplotypes/${sample}.hap1.fa contigsFasta
cat /home/haplotypes/${sample}.hap2.fa >> contigsFasta

cp  ../sample-hap1/yahs.out_scaffolds_final.agp scaffold.joint.agp
cat ../sample-hap2/yahs.out_scaffolds_final.agp >> scaffold.joint.agp

#Align Hi-C data to the assembly, remove PCR duplicates and filter out secondary and supplementary alignments
bwa index contigsFasta
bwa mem -5SP contigsFasta r1Reads r2Reads -t 9 | samblaster | samtools view - -@ 9 -S 
-h -b -F 3340 -o HiC.bam

rm -rf HiC.filtered.bam
#Filter the alignments with MAPQ 0 (mapping quality =1) and NM 3 (edit distance < 3)
/home/software/HapHiC/utils/filter_bam HiC.bam 1 --nm 3 --threads 9 |samtools view - -b -@ 9 -o HiC.filtered.bam

samtools faidx contigsFasta
/home/software/HapHiC/scripts/../utils/juicer pre -a -q 1 -o out_JBAT HiC.filtered.bam scaffold.joint.agp contigsFasta.fai >out_JBAT.log 2>&1
```

### Manual scaffold adjustment

The Hi-C contact matrix was manually refined using Juicebox. Adjustments mainly involved correcting contig orientations and orders. Contigs with phasing errors or mis-assemblies were interrupted only when necessary, and a few abnormal signal regions were removed to minimize potential errors.

## Assembly Result Visualization

```shell
/public/home/user/software/HapHiC/haphic plot scaffold.joint.agp HiC.bam
```

