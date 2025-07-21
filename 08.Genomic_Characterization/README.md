# Genomic characterization of haplotypes

- [Genome-wide heterozygosity](#genome-wide-heterozygosity)
- [Haplotype SNP/SV-PCA](#haplotype-snpsv-pca)
- [Haplotype combination](#haplotype-combination)

------

## Genome-wide heterozygosity

To profile heterozygous regions within four wild accessions and 54 cassava cultivars, we implemented a 200-kb non-overlapping sliding window approach to classify genomic segments as either homozygous or heterozygous. Heterozygous variants were identified by alignment of HiFi reads to haplotype1 for each accession. Genomic windows exhibiting a heterozygous SNP density exceeding 2 SNPs/kb were classified as heterozygous segments; all other windows were designated as homozygous. The total heterozygous genome size for each individual was then quantified by summing the lengths of all identified heterozygous windows.

```shell
vcftools --vcf ../all.miss0.7.biallic.vcf --recode --indv ${sample} --stdout > ${sample}.vcf
grep -v '0|0' ${sample}.vcf | grep -v '1|1' | grep -v './.' > ${sample}_het.vcf
vcftools --vcf ${sample}_het.vcf --SNPdensity 200000 --out ${sample}
awk -v OFS='\t' '$3 = $2 + 200000' ${sample}.snpden > ${sample}.snpden.2
python merge_block.py  -i ${sample}.snpden.2 -t 2 > ${sample}.2.het_hom_region.txt #The merge_block.py script is included in the same directory as README.md
```

------

## Haplotype SNP/SV-PCA

To dissect haplotype combinatorial patterns in cultivated cassava, we performed principal component analysis (PCA) on paired haplotypes from 54 accessions. We conducted independent PCA analyses for each of the 18 chromosomes, employing both SNP and SV datasets. 

```shell
vcftools --vcf variant.vcf --keep all_CUL_sample.txt --recode --recode-INFO-all --stdout > CUL_variant.vcf

# PCA
for i in `cat chr`  #chr: A file containing all chromosome names in the reference genome.
do
bcftools view -h CUL_variant.vcf > ${i}.vcf
cat CUL_variant.vcf | grep -v '#' | grep "${i}"  > ${i}_tmp.vcf
cat ${i}_tmp.vcf >> ${i}.vcf
plink --vcf ${i}.vcf --allow-extra-chr  --double-id --make-bed --out ${i}
plink --bfile ${i} --pca 10 --out ${i}_results
echo -e "${i}_PCA1\t${i}_PCA2" > ${i}.PCA1.PCA2.txt
awk '{print $3 "\t" $4}' ${i}_results.eigenvec >> ${i}.PCA1.PCA2.txt
done

# Statistics
for  i in `cat chr`
do
echo -e "${i}_PCA1\t${i}_PCA2" > ${i}.PCA1.PCA2.txt
awk '{print $3 "\t" $4}' ${i}_results.eigenvec >> ${i}.PCA1.PCA2.txt
echo "${i}" >> contrubution.txt
head -n 2 chr18_results.eigenval >> contrubution.txt
done

awk '{print $1}' chr17_results.eigenvec > ample.PCA1.PCA2.txt
paste *.PCA1.PCA2.txt > all.PCA1.PCA2.txt
```



------

## Haplotype combination

We constructed a chromosome-scale haplotype sharing landscape for the cassava pan-genome based on sequence similarity but without considering recombination events. Our methodology comprised the following steps: Initially, we computed pairwise haplotype sequence divergence across the genome using a 100-kb non-overlapping sliding window applied to phased single nucleotide polymorphism (SNP) datasets. Subsequently, we generated inter-accession genetic distance matrices employing the scipy.spatial.distance module and performed hierarchical clustering via the scipy.cluster.hierarchy module. The optimal number of discrete haplotype clusters was then determined by maximizing the silhouette coefficient, with a predefined upper bound of 10 clusters to mitigate potential overfitting. Following the integration of genome-wide window data, the resulting haplotype clusters were ranked and visually differentiated by color based on their constituent haplotype counts. The two haplotypes from each individual are interconnected by a dashed line to visualize their combination.

```shell
# Step 1: Data Preprocessing Pipeline
# Split VCF into 100kb windows
python 01_split_VCF_into_100kb_window.py -i CUL_variant.vcf -o prefix_dir

cd prefix_dir

# Process each window sequentially
for i in *.vcf; do
    # Convert VCF to TSV format
    python 02_transform_VCF_to_TSV.py -i $i -o ${i}.tsv
    
    # Perform initial clustering
    python 03_cluster_based_on_TSV.py -i ${i}.tsv -o ${i}.cluster
    
    # Transform clustering orientation from columns to rows
    python 04_clustering_by_columns_to_by_rows.py -i ${i}.cluster -o ${i}.rows
done

# Step 2: Result Consolidation and Visualization
# Merge all cluster results chromosome-wise by genomic position
# Then generate visualization plots

```

