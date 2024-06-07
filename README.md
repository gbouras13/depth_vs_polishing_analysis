# How Low Can You Go? Short-read polishing of Oxford Nanopore bacterial genome assemblies Code Repository

This repository holds code and data for this manuscript:<br>
[Bouras G, Judd LM, Edwards RA, Vreugde S, Stinear TP, Wick RR. How low can you go? Short-read polishing of Oxford Nanopore bacterial genome assemblies. Microbial Genomics. 2024. doi:10.1099/mgen.0.001254.](https://doi.org/10.1099/mgen.0.001254)



Contents:
* `figures`: contains all of the manuscript's main and supplementary figures along with their captions.
* `supp_tables.xlsx`: contains the paper's supplementary tables.
* `ont_assemblies`: contains the ONT-only [Trycycler](https://github.com/rrwick/Trycycler) assemblies used as a starting point for polishing.
* `reference_assemblies`: contains the polished and manually curated assemblies used as a ground truth.
* `pypolca_example_plot`: contains code to simulate reads, errors and make (Figure 1).
* `main_analysis`: contains the read subsampling, polishing and plotting commands for the main analysis (Figures 2, S2, S3 and S7).
* `errors_in_repeats`: contains the details of the errors-in-repeats analysis (Figure S1).
* `long_homopolymer`: contains the details of the long-homopolymer analysis (Figure S4).
* `error_characterisation`: contains the detailed error characterisation of the 37 existing errors and all polisher introduced errors (Table S2 and Figures S5 and S6).
* `hybracter_analysis`: contains the read subsampling assembly and plotting commands for the [Hybracter](https://github.com/gbouras13/hybracter) analysis in the paper (Figures S8 and S9).
* `reference_chromosome_assemblies_hybracter`: contains the polished and manually curated assemblies used as a ground truth, chromosomes only. Used for the Hybracter analysis.
* `low_quality_draft`: contains the details of the polishing analysis using low-quality draft assemblies (Figure S10).
* `parameter_sweep`: contains the details of the low-depth parameter sweep analysis (Table S6).
* `compare_assemblies.py`: assembly comparison script used for counting/characterising errors.


ONT and Illumina reads are not included in this repository due to size, but they can be found on SRA:

| Genome                                  | ONT reads                                                         | Illumina reads                                                    |
|-----------------------------------------|-------------------------------------------------------------------|-------------------------------------------------------------------|
| _Campylobacter jejuni_ (ATCC-33560)     | [SRR27638397](https://www.ncbi.nlm.nih.gov/sra/?term=SRR27638397) | [SRR26899120](https://www.ncbi.nlm.nih.gov/sra/?term=SRR26899120) |
| _Campylobacter lari_ (ATCC-35221)       | [SRR27638396](https://www.ncbi.nlm.nih.gov/sra/?term=SRR27638396) | [SRR26899115](https://www.ncbi.nlm.nih.gov/sra/?term=SRR26899115) |
| _Escherichia coli_ (ATCC-25922)         | [SRR27638398](https://www.ncbi.nlm.nih.gov/sra/?term=SRR27638398) | [SRR26899128](https://www.ncbi.nlm.nih.gov/sra/?term=SRR26899128) |
| _Listeria ivanovii_ (ATCC-19119)        | [SRR27638399](https://www.ncbi.nlm.nih.gov/sra/?term=SRR27638399) | [SRR26899136](https://www.ncbi.nlm.nih.gov/sra/?term=SRR26899136) |
| _Listeria monocytogenes_ (ATCC-BAA-679) | [SRR27638394](https://www.ncbi.nlm.nih.gov/sra/?term=SRR27638394) | [SRR26899101](https://www.ncbi.nlm.nih.gov/sra/?term=SRR26899101) |
| _Listeria welshimeri_ (ATCC-35897)      | [SRR27638395](https://www.ncbi.nlm.nih.gov/sra/?term=SRR27638395) | [SRR26899109](https://www.ncbi.nlm.nih.gov/sra/?term=SRR26899109) |
| _Salmonella enterica_ (ATCC-10708)      | [SRR27638402](https://www.ncbi.nlm.nih.gov/sra/?term=SRR27638402) | [SRR26899135](https://www.ncbi.nlm.nih.gov/sra/?term=SRR26899135) |
| _Vibrio cholerae_ (ATCC-14035)          | [SRR27638401](https://www.ncbi.nlm.nih.gov/sra/?term=SRR27638401) | [SRR26899095](https://www.ncbi.nlm.nih.gov/sra/?term=SRR26899095) |
| _Vibrio parahaemolyticus_ (ATCC-17802)  | [SRR27638400](https://www.ncbi.nlm.nih.gov/sra/?term=SRR27638400) | [SRR26899141](https://www.ncbi.nlm.nih.gov/sra/?term=SRR26899141) |

These are easily downloaded using the [`fastq-dl`](https://github.com/rpetit3/fastq-dl) program e.g.

```
CPUS=16
fastq-dl --accession SRR27638397 --cpus $CPUS
fastq-dl --accession SRR27638396 --cpus $CPUS
fastq-dl --accession SRR27638398 --cpus $CPUS
fastq-dl --accession SRR27638399 --cpus $CPUS
fastq-dl --accession SRR27638394 --cpus $CPUS
fastq-dl --accession SRR27638395 --cpus $CPUS
fastq-dl --accession SRR27638402 --cpus $CPUS
fastq-dl --accession SRR27638401 --cpus $CPUS
fastq-dl --accession SRR27638400 --cpus $CPUS
```
