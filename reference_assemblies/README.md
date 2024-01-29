This directory contains references assemblies for each of the genomes used in this study. They assembled using Trycycler followed by Illumina polishing (using various tools) and manual curation. They should be completely error-free.

The _Escherichia coli_ genome contains some structural heterogeneity in its largest plasmid - there is a 3 kbp sequence which occurs in both orientations:
```
Structure 1:
========>>>>>>>>>>========

Structure 2:
========<<<<<<<<<<========
```
So there are two valid ways to assemble this plasmid. The reads (`ATCC_25922_Escherichia_coli_reads.fasta.gz`) support both orientations in a near 50:50 split, so assemblers may choose either option. `ATCC_25922_Escherichia_coli.fasta.gz` contains the main reference assembly and `ATCC_25922_Escherichia_coli_alt.fasta.gz` contains the other structural variant in that plasmid.

Similarly, the Vibrio cholerae_ genome contains some structural heterogeneity in its largest chromosome - there is a 34 kbp inverted repeat with 44 kbp of sequence in the middle, and the middle sequence can occur in both possible orientations:
```
Structure 1:
========-------->>>>>>>>>>--------========
         repeat            repeat

Structure 2:
========--------<<<<<<<<<<--------========
         repeat            repeat
```
So there are two valid ways to assemble this chromosome. Since the repeat is quite long, only 10 ONT reads directly inform this structure (`ATCC_14035_Vibrio_cholerae_reads.fasta.gz`). Eight reads support structure 1, and two reads support structure 2, but these two reads (`401525da` and `6d051173`) appear to be a duplex pair, so they are effectively only one read. However, this is a very long read (over 80 kbp) and so assemblers may give it priority when constructing the assembled sequence, and this sometimes leads to misassemblies. `ATCC_14035_Vibrio_cholerae.fasta.gz` contains the main reference assembly (structure 1) and `ATCC_14035_Vibrio_cholerae_alt.fasta.gz` contains the other variant (structure 2).
