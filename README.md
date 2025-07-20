<h1 align="center">Cassava Pangenome Project</h1>

![cassava_pangenome](./cassava_pangenome.png)

## Project Introduction

This repository hosts the code and resources for constructing a pangenome graph of cassava (*Manihot esculenta*), a clonally propagated staple crop feeding millions worldwide. By analyzing **115 diverse haplotypes**, we built a pangenome **2.6× larger** than the linear reference, revealing how clonal propagation preserves heterozygosity and structural variation (SV) through unique evolutionary mechanisms. 

Explore the data, tools, and pipelines to study SV, haplotype evolution, and clonal genome dynamics. 

------

## Document

### [1. Haplotype-Resolved Assembly](01.Haplotype_Resolved_Assembly/README.md)
- Reference genome assembly
- Haplotype-resolved assembly
- Contig sequence filtering
- Haplotype-resolved scaffolding
- Assembly visualization

### [2. Genome Quality Evaluation](02.Genome_Quality_Evaluation/README.md)
- Contiguity metrics (N50/L50)
- Completeness assessment
  - BUSCO analysis
  - Merqury evaluation
- Assembly correctness
  - Switch/Harming errors
  - Mis-assembly detection

### [3. Genome Annotation](03.Genome_Annotation/README.md)
- Transposable elements
- Gene prediction

### [4. Sequence Divergence](04.Sequence_Divergence/README.md)
- Mash distance
- Mapping rates

### [5. Graph Pangenome Construction](05.Graph_Pangenome_Construction/README.md)
- Graph construction pipeline
  - Wfmash alignment
  - Seqwish graph
  - PGGB processing
- Graph analysis
  - Growth curve
  - Non-reference nodes

### [6. Structural Variant Identification](06.Structural_Variants/README.md)
- SV detection
  - Haplotype alignment
  - SVIM-asm pipeline
  - SYRI inversions
  - Jasmine merging
- Multiallelic variants

### [7. Breeding History Analysis](07.Breeding_History/README.md)
- Kinship analysis
- IBD analysis
- Effective population size

### [8. Haplotype Characterization](08.Haplotype_Characterization/README.md)
- Heterozygosity analysis
- Haplotype PCA
- Haplotype combinations

### [9. Genetic Burden Analysis](09.Genetic_Burden/README.md)
- Minor-frequency SVs
- S1 progeny simulation
  - Selfing simulations
  - Hybrid simulations

### [10. Gene Duplication](10.Gene_Duplication/README.md)
- Orthogroups analysis
- Duplication events

### [11. Paleotetraploidy Analysis](11.Paleotetraploid/README.md)
- WGD detection
- Karyotype reconstruction
- Core genes identification

### [12. Biased Fractionation](12.Biased_Fractionation/README.md)
- Gene retention analysis
- Unbalanced regions
- LF/HF comparisons

### [13. Translocation Identification](13.Translocation/README.md)
- Network analysis
- Hi-C validation
- Breakpoint analysis
- Synteny analysis

### [14. Translocation Characterization](14.Translocation_Characterization/README.md)
- Timing estimation
- Centromere analysis
- Diversity analysis
- Linkage analysis

------

## Citation

If you use this toolkit in your research, please cite:

```
XXX, XXX., et al. (2025). XXX XXX XXX XXX, 15(3), 234-245.
```

------

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
