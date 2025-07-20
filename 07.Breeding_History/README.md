# Cassava breeding history 

- [Kinship analysis](#kinship-analysis)
- [Identity-by-descent (IBD) analysis](#identity-by-descent-ibd-analysis)
- [Effective population size (Ne)](#effective-population-size-ne)

------

## Kinship analysis

We performed kinship analysis using the [KING](https://www.chen.kingrelatedness.com/) software (v2.3.2) with the --kinship option based on the phased SNP map of cassava. Sample relationships were classified according to the InfType results: duplicate/monozygotic (Dup/MZ) samples were merged; parent-offspring (PO) and full sibling (FS) relationships were grouped as first-degree relatives. 

```

```



------

## Identity-by-descent (IBD) analysis

Our assembly data achieved chromosome-level phasing. We aligned haplotype-resolved genomes of 54 cultivated cassava individuals (108 haplotypes) to the reference genome AM560, generating a phased SNP map. Variant sites were filtered to retain only those with a mapping rate >90%. We then processed the filtered variants using [BEAGLE](https://github.com/beagle-dev/beagle-lib) (v27Feb25.75f) with default parameters (assuming a recombination rate of 1 cM/Mb) to: (1) impute missing genotypes and (2) refine haplotype phasing to reduce errors. Subsequently, we inferred IBD segments using [hap-ibd.jar](https://hpc.nih.gov/apps/Hap-IBD.html) (v1.0.0) with this phased SNP dataset with default parameters. As recommended, we focused on larger IBD fragments because they suggest recent hybridization events, and thus only retained IBD fragment pairs larger than 3 Mb. We then quantified the number of IBD segments and their cumulative length for each haplotype pair. 

```

```



------

## Effective population size (*Ne*)

We reconstruct recent effective population size dynamics (past ~100 generations) based on the length distribution of IBD segments shared among individuals. IBD segment length is inversely correlated with the time to the most recent common ancestor (TMRCA); shorter segments reflect older coalescence events, while longer segments indicate more recent shared ancestry. Specifically, we processed the data using [IBDNe](https://faculty.washington.edu/browning/ibdne.html) (v23Apr20.ae9), and retained IBD segments >1 Mb in length. The analysis incorporated an assuming recombination rate parameter of 1 cM/Mb. We then generated reconstructions of demographic history, producing temporal *Ne* estimates with 95% confidence intervals spanning the past 100 generations of evolutionary history.

```

```

