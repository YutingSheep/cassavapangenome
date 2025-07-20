# Genetic Burden

- [Minor-frequency SVs](#minor-frequency-svs)
- [Simulated S1 progeny](#simulated-s1-progeny)
  - [The input VCF format](#the-input-vcf-format-for-the-simulated-population-is-as-follows)
  - [Selfing simulations](#selfing-simulations)
  - [Hybrid simulations](#hybrid-simulations)

------

## Minor-frequency SVs

------

We focused our investigation on derived SVs exhibiting minor allele frequencies (MAF) below 0.05 in the population. This threshold was selected to capture evolutionarily recent or potentially deleterious variants that might be under selection pressure. From our dataset of 108 cultivated haplotypes, we established a count threshold of 10 (approximately 9.3% of samples) to define minor-frequency SVs. We performed genotype analysis across all haplotypes, categorizing them as either heterozygous or homozygous derived states. We then constructed distributions for: (1) all derived minor-frequency SVs, (2) heterozygous SVs, and (3) homozygous derived SVs. 

```shell
mkdir result
cd result
input_list=($(cat ../input_check.txt))

process_sample() {
    local i="$1"
    local tmp_result=$(mktemp)

    grep "_${i}_" ../all_sample.txt > "${i}.txt"
    
    vcftools --vcf ../output_filtered.vcf \
             --keep "${i}.txt" \
             --recode --recode-INFO-all \
             --stdout > "${i}.vcf" 2>/dev/null
    
    awk '{print $1 "\t" $2 "\t" $10 "\t" $11}' "${i}.vcf" | 
        grep '0/0' | grep '1/1' > "${i}_het.vcf"
        awk '{print $1 "\t" $2 "\t" $10 "\t" $11}' "${i}.vcf" | grep -v '0/0' > "${i}_hom.vcf"

    echo "${i}" >> "$tmp_result"
    while read k; do
        grep_count1=$(grep "${k}" "${i}_het.vcf" | awk '{print $3}' | grep '1/1' | wc -l)
        grep_count2=$(grep "${k}" "${i}_het.vcf" | awk '{print $4}' | grep '1/1' | wc -l)
        echo "$grep_count1" >> "$tmp_result"
        echo "$grep_count2" >> "$tmp_result"
    done < ../list
    
    flock all_result.txt -c "cat '$tmp_result' >> all_result.txt"
    #rm "$tmp_result" "${i}.txt" "${i}.vcf" "${i}_het.vcf"
}

export -f process_sample
parallel -j 8 process_sample ::: "${input_list[@]}"
cd ..
```



## Simulated S1 progeny

To investigate the distribution of derived SVs in selfed (or hybrid) S1 progeny, we leveraged genomic data from 54 cultivated accessions. Under the idealized assumption of no intra-chromosomal recombination, we modeled the S1 generation as a random combination of the two haplotypes across 18 chromosomes. This approach allowed us to simulate the inheritance of derived SVs while maintaining linkage blocks intact. 

### The input VCF format for the simulated population is as follows:

| #CHROM | POS  | 0_BRA001554_hap1 | 1_BRA001554_hap2 | 2_BRA117315_hap1 | 3_BRA117315_hap2 | 4_s100_hap1 | 5_s100_hap2 |
| :----: | :--: | :--------------: | :--------------: | :--------------: | :--------------: | :---------: | ----------- |
| chr01  | 158  |       0/0        |       0/0        |       0/0        |       0/0        |     0/0     | 0/0         |
| chr01  | 326  |       0/0        |       0/0        |       0/0        |       0/0        |     0/0     | 0/0         |
| chr01  | 341  |       0/0        |       0/0        |       0/0        |       0/0        |     0/0     | 0/0         |
| chr01  | 434  |       0/0        |       0/0        |       0/0        |       0/0        |     0/0     | 0/0         |
| chr01  | 510  |       0/0        |       0/0        |       0/0        |       0/0        |     0/0     | 0/0         |
| chr01  | 528  |       0/0        |       0/0        |       1/1        |       1/1        |     0/0     | 0/0         |
| chr01  | 807  |       0/0        |       0/0        |       0/0        |       0/0        |     0/0     | 0/0         |
| chr01  | 825  |       0/0        |       0/0        |       0/0        |       0/0        |     0/0     | 0/0         |
| chr01  | 873  |       0/0        |       0/0        |       0/0        |       0/0        |     0/0     | 0/0         |

### Selfing simulations

For selfing simulations, we computed the expected distribution of derived SVs in S1 offspring by randomly pairing haplotypes from the same individual, thereby simulating the consequences of inbreeding. 

```shell
##The simulation_selfing.py script is included in the same directory as README.md
mkdir -p ./result/ && \
seq 1 2000 | parallel -j 30 'python simulation_selfing.py -i "unfolded_2hap_DEL_INS_for_simulation.vcf" -o "./result/{}_sim.txt"'
```

### Hybrid simulations

Similarly, for hybrid simulations, we modeled crosses between distinct parental haplotypes to estimate the distribution of derived SVs in hybrid progeny. 

```shell
##The simulation_hybrid.py script is included in the same directory as README.md
mkdir -p ./result/ && \
seq 1 2000 | parallel -j 30 'python simulation_hybrid.py -i "unfolded_2hap_DEL_INS_for_simulation.vcf" -o "./result/{}_sim.txt"'
```

