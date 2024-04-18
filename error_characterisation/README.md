This directory holds code and analysis for the error characterisation analysis.

There are 37 total errors between the ground truth reference genomes and the long-read only Trycycler assemblies. 

One error is characterised in detail in the manuscript (Supplementary Figure 1) - the extra-long homopolymer in ATCC 10708 _Salmonella enterica_ , the deep dive of which can be found in the `long_homopolymer` repository.

The exact loci and representation of all errors can be found in the `ont_assemblies` directory in the `.errors` files - the top row represents the long-read only Trycycler assembly, while the bottom row is the reference genome.

# Characterisation

A characterisation of each of these errors can be found in `all_errors.xlsx`.

There are 32 error loci containing 37 total errors - some deletion have multiple errors (see column D `Indel Size (bp)` for those of size larger than 1).

Of the 32 error loci, 17 are SNPs, 14 are deletions (all in homopolymers) and 1 is an insertion.

13/17 SNPs occur in ATCC 35221 _Campylobacter lari_ and are all C->T or G->A changes in the GATC motif.

### Methylation

Therefore, we decided to conduct a methylation analysis of _C. lari_ with [MicrobeMod](https://github.com/cultivarium/MicrobeMod). See the `methylation` directory for more details. As can be seen, there was at least one 6mA methylation near every SNP, with most occuring within 2 bases.

Therefore, it appears possible that long-only SNP errors in this genome are due to nearby methylation. We did not do further analysis here 

How much does MicrobeMod rely on correct read basecalls at its motifs? I.e if the Dorado is regularly screwing up at some sites (which I think it must be if we're getting assembly errors), does that make Dorado's methylation calls also screwed up? And then does that impact MicrobeMod's ability to correctly call a motif?
It might be worth rerunning MicroMod with a lower value for --percent_methylation_cutoff. I can see that some of the modified bases near our errors are quite close to the default value of 0.66. Perhaps this is related to the above point - bad basecalls might cause MicrobeMod to think it's a lower percentage modified than it actually is.
GATC is its own reverse complement, so I'd expect to see both strands methylated at GATC sites. This is the case for some but not all of our errors. I wonder if a lower percent methylation threshold would make it true for all of them?
We could test the above points by filtering out any reads which have errors near our assembly errors (I'm assuming there are plenty of error-free reads at these sites). My hypothesis is that MicrobeMod would then report a very high methylation percent for both strands, and the picture would become more clear.

# Error Analysis By Depth

We also looked at the ability of each polishing tool to fix each of the 37 errors vs whether they introduced new errors, at higher (25x+) depths. We didn't do depths lower than 25x as it is obvious that the polishers other than Pypolca-careful and Polypolish in both modes commonly introduce hundreds or thousands of errors.

For example, we wanted to know that if a polisher left 10 remaining errors, was this because it fixed 27/37 original errors without introducing any, or did it fix all 37 and introduce an extra 10?

To do this, we used the `compare_assemblies.py` output (explained in the `main_analysis` directory) for every assembly from 25.0x to 50.0x depth. We then used the script `analyse_errors.py` to determine whether the remaining errors were part of the original 37, or introduced, based on the error region.

```
python analyse_errors.py -i assemblies -o analysed_output  --min_depth 25.0  --max_depth 50.0 -c ont_error_list.csv -f
```

Specifically the `assemblies` directory looked like 

```
├── assemblies
│   ├── 25.0
│   │   ├── ATCC_10708_Salmonella_enterica
│   │   │   ├── fmlrc2.errors
│   │   │   ├── hypo.errors
│   │   │   ├── nextpolish.errors
│   │   │   ├── pilon.errors
│   │   │   ├── polypolish-careful.errors
│   │   │   ├── polypolish.errors
│   │   │   ├── pypolca-careful.errors
│   │   │   ├── pypolca.errors
│   │   ├── ATCC_14035_Vibrio_cholerae
│   │   │   ├── fmlrc2.errors
│   │   │   ├── hypo.errors
│   │   │   ├── nextpolish.errors
│   │   │   ├── pilon.errors
│   │   │   ├── polypolish-careful.errors
│   │   │   ├── polypolish.errors
│   │   │   ├── pypolca-careful.errors
│   │   │   ├── pypolca.errors
│   │   ├── ...
│   ├── 25.1
│   │   ├── ...
│   ├── 25.2
│   ├── ...

```

## Per Error

An analysis of every error region identified by `compare_assemblies.py` for every genome and depth can be found in `analysed_output`/`per_error_output.csv`. The `existing` column is True if it is one of the 37 original errors. Otherwise, that column if False and the error is introduced. `error_count` is more than 1 if there were multiple errors in the genomic region:

 This error region would have 3 errors (indicated by an *) - note: the top row represents the long-read only Trycycler assembly, while the bottom row is the reference genome.

```
chromosome 1938765-1938801: GTGCAGGCTGGTCACAGTTGACTCGATTTGCGTCATC
chromosome 1938766-1938802: GTGCAGGCTGGTCACGGTCGATTCGATTTGCGTCATC
                                           *  *  *       
```

## Aggregated

`analysed_output`/`aggregated_output.csv` contains the sum of each error type per genome and depth: SNP, Insertion, Deletion or Mixed (more than one of the above in the same genome region), along with sums for existing and introduced errors.

**Note: the top row represents the long-read only Trycycler assembly, while the bottom row is the reference genome.**

Examples of each type are as follows:

### SNP

```
chromosome 1938765-1938801: GTGCAGGCTGGTCACAGTTGACTCGATTTGCGTCATC
chromosome 1938766-1938802: GTGCAGGCTGGTCACGGTCGATTCGATTTGCGTCATC
                                           *  *  *       
```

### Deletion

```

chromosome 548169-548198: TTTTTTTTTTTTTTT---GAAAGTATGCAACAA
chromosome 548170-548202: TTTTTTTTTTTTTTTTTTGAAAGTATGCAACAA
                                        ***      

```

### Insertion

```
chromosome 550320-550350: CGGATCCCCCCCCCCCGATAGCGTGAGCTGA
chromosome 550320-550349: CGGATCCCCCCCCCC-GATAGCGTGAGCTGA
                                        *              
```

### Mixed 

* This has SNPs, deletions and an insertion (the last error).
* Therefore this is denoted as a 'Mixed' error region - these are the worst kind introduced by polishing tools

```
chromosome 3119155-3119203: CCCCCAGCCTAGCTGTAA-TGC-CAGT-CAGTTAA-GCAGCCAGTCAGTTGCT
chromosome 3119156-3119207: CCCCCAGCCTAGCTGGGGGTTTTCTGTGCACAAAAAG-AGCCAGTCAGTTGCT
                                           **** *** *  *  ***  * *               
```

## Conclusions

* Polypolish (both careful and default), pypolca-careful rarely introduce errors (we knew this from the main analysis!)
    * Polypolish-careful introduced errors twice, Polypolish-default eight times and Pypolca-careful 54 times across all depths
* Interestingly, the other 4 polishers could in fact polish almost all existing errors (the min existing errors remaining were 0 or 1 for all 8 polishers) but they commonly introduced errors
* Nextpolish and especially HyPo consistently introduced errors even up to 50x
* With fmlrc2, there was a clear relationship with depth - increasing depth led to fewer introduced errors
* Pilon had extremely high variance - at some depths it performed well, at others it introduced hundreds of errors

To recreate Supplementary Figures 8 and 9, please run `plots.Rmd`
