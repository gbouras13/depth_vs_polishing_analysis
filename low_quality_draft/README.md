This analysis repeats the main analysis, but instead of using high quality draft assemblies (Trycycler from sup reads), it uses lower quality draft assemblies (Trycycler from fast reads).



## Starting files

For each genome, I'm starting with the draft (ONT-only) and reference (manually curated ground truth):
```
drafts_low_quality/ATCC_10708_Salmonella_enterica.fasta
drafts_low_quality/ATCC_14035_Vibrio_cholerae.fasta
drafts_low_quality/ATCC_17802_Vibrio_parahaemolyticus.fasta
drafts_low_quality/ATCC_19119_Listeria_ivanovii.fasta
drafts_low_quality/ATCC_25922_Escherichia_coli.fasta
drafts_low_quality/ATCC_33560_Campylobacter_jejuni.fasta
drafts_low_quality/ATCC_35221_Campylobacter_lari.fasta
drafts_low_quality/ATCC_35897_Listeria_welshimeri.fasta
drafts_low_quality/ATCC_BAA-679_Listeria_monocytogenes.fasta

references/ATCC_10708_Salmonella_enterica.fasta
references/ATCC_14035_Vibrio_cholerae.fasta
references/ATCC_17802_Vibrio_parahaemolyticus.fasta
references/ATCC_19119_Listeria_ivanovii.fasta
references/ATCC_25922_Escherichia_coli.fasta
references/ATCC_33560_Campylobacter_jejuni.fasta
references/ATCC_35221_Campylobacter_lari.fasta
references/ATCC_35897_Listeria_welshimeri.fasta
references/ATCC_BAA-679_Listeria_monocytogenes.fasta
```



## Set variables

The directories in `base_dir`, `read_dir` and `nextpolish_dir` are for my system - change as appropriate for your system. Genome sizes were gotten from the ground truth reference genomes.

```bash
base_dir=/home/wickr/2024-01_low-depth_Illumina_polishing
read_dir=/home/wickr/2023-10_PulseNet_ATCC_assemblies
nextpolish_dir=/home/wickr/programs/NextPolish

genomes=(
    "ATCC_10708_Salmonella_enterica"
    "ATCC_14035_Vibrio_cholerae"
    "ATCC_17802_Vibrio_parahaemolyticus"
    "ATCC_19119_Listeria_ivanovii"
    "ATCC_25922_Escherichia_coli"
    "ATCC_33560_Campylobacter_jejuni"
    "ATCC_35221_Campylobacter_lari"
    "ATCC_35897_Listeria_welshimeri"
    "ATCC_BAA-679_Listeria_monocytogenes"
)

declare -A genome_sizes
genome_sizes["ATCC_10708_Salmonella_enterica"]=4801704
genome_sizes["ATCC_14035_Vibrio_cholerae"]=4142374
genome_sizes["ATCC_17802_Vibrio_parahaemolyticus"]=5147091
genome_sizes["ATCC_19119_Listeria_ivanovii"]=2919549
genome_sizes["ATCC_25922_Escherichia_coli"]=5204893
genome_sizes["ATCC_33560_Campylobacter_jejuni"]=1768448
genome_sizes["ATCC_35221_Campylobacter_lari"]=1513368
genome_sizes["ATCC_35897_Listeria_welshimeri"]=2814137
genome_sizes["ATCC_BAA-679_Listeria_monocytogenes"]=2944530
```



## Build indices

Build BWA and samtools indices for the draft genomes:
```bash
cd "$base_dir"/drafts_low_quality
for a in "${genomes[@]}"; do
    bwa index "$a".fasta
    samtools faidx "$a".fasta
done
```



## Error counts in draft genomes

```bash
cd "$base_dir"/drafts_low_quality

for a in "${genomes[@]}"; do
    ref="$base_dir"/references/"$a".fasta
    compare_assemblies.py --aligner edlib "$ref" "$a".fasta > "$a".errors
done

for a in "${genomes[@]}"; do
    printf "$a\t"
    cat "$a".errors | grep -o "*" | wc -l
done
```

```
ATCC_10708_Salmonella_enterica:      4132
ATCC_14035_Vibrio_cholerae:          4664
ATCC_17802_Vibrio_parahaemolyticus:  4202
ATCC_19119_Listeria_ivanovii:         936
ATCC_25922_Escherichia_coli:         5510
ATCC_33560_Campylobacter_jejuni:     2538
ATCC_35221_Campylobacter_lari:       2739
ATCC_35897_Listeria_welshimeri:       668
ATCC_BAA-679_Listeria_monocytogenes:  795
-----------------------------------------
total                               26185
```

This is roughly 3 orders of magnitude more errors than in the main analysis (Q30 vs Q60).



## Run individual polishers

Polypolish v0.6.0 (both defaults and --careful):
```bash
conda activate assembly

for d in $(seq -f "%04.1f" 0.1 0.1 50); do
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        draft="$base_dir"/drafts_low_quality/"$a".fasta

        bwa mem -t 24 -a "$draft" reads_1.fastq.gz > alignments_1.sam 2>> bwa.txt
        bwa mem -t 24 -a "$draft" reads_2.fastq.gz > alignments_2.sam 2>> bwa.txt
        polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam 2> polypolish_filter.txt
        polypolish polish "$draft" filtered_1.sam filtered_2.sam 1> polypolish.fasta 2> polypolish.txt
        polypolish polish --careful "$draft" filtered_1.sam filtered_2.sam 1> polypolish-careful.fasta 2> polypolish-careful.txt
        rm *.sam

        gzip polypolish.fasta polypolish-careful.fasta
    done
done
```

pypolca v0.3.0 (both defaults and --careful):
```bash
conda activate assembly

for d in $(seq -f "%04.1f" 0.1 0.1 50); do
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        draft="$base_dir"/drafts_low_quality/"$a".fasta

        pypolca run -a "$draft" -1 reads_1.fastq.gz -2 reads_2.fastq.gz -t 24 -o pypolca > pypolca.txt 2>&1
        seqtk seq pypolca/pypolca_corrected.fasta > pypolca.fasta
        rm -r pypolca

        pypolca run --careful -a "$draft" -1 reads_1.fastq.gz -2 reads_2.fastq.gz -t 24 -o pypolca > pypolca-careful.txt 2>&1
        seqtk seq pypolca/pypolca_corrected.fasta > pypolca-careful.fasta
        rm -r pypolca

        gzip pypolca.fasta pypolca-careful.fasta
    done
done
```

HyPo v1.0.3:
```bash
conda activate hypo

for d in $(seq -f "%04.1f" 0.1 0.1 50); do 
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        draft="$base_dir"/drafts_low_quality/"$a".fasta
        genome_size=${genome_sizes["$a"]}
        integer_depth=$(printf "%.0f" "$d"); if [[ $integer_depth == 0 ]]; then integer_depth=1; fi

        bwa mem -t 24 "$draft" reads_1.fastq.gz reads_2.fastq.gz 2> /dev/null | samtools sort > hypo_alignments.bam; samtools index hypo_alignments.bam
        echo -e "reads_1.fastq.gz\nreads_2.fastq.gz" > read_filenames.txt
        while [[ ! -f "hypo.fasta" ]]; do  # HyPo randomly fails a lot when depth <5, so repeated tries are necessary.
            hypo -d "$draft" -r @read_filenames.txt -s "$genome_size" -c "$integer_depth" -b hypo_alignments.bam -t 24 -o hypo.fasta > hypo.txt 2>&1
            rm -r aux
        done
        rm hypo_alignments.bam hypo_alignments.bam.bai read_filenames.txt

        gzip hypo.fasta
    done
done
```

FMLRC2 v0.1.8:
```bash
conda activate assembly

for d in $(seq -f "%04.1f" 0.1 0.1 50); do 
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        draft="$base_dir"/drafts_low_quality/"$a".fasta

        gunzip -c reads_1.fastq.gz reads_2.fastq.gz | awk 'NR % 4 == 2' | tr NT TN | ropebwt2 -LR | tr NT TN | fmlrc2-convert comp_msbwt.npy 2> fmlrc2.txt
        fmlrc2 -t 24 comp_msbwt.npy "$draft" fmlrc2.fasta 2>> fmlrc2.txt
        rm comp_msbwt.npy

        gzip fmlrc2.fasta
    done
done
```

NextPolish v1.4.1 (1](https://nextpolish.readthedocs.io/en/latest/TUTORIAL.html)):
```bash
conda activate nextdenovo

for d in $(seq -f "%04.1f" 0.1 0.1 50); do 
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        draft="$base_dir"/drafts_low_quality/"$a".fasta

        bwa mem -t 24 "$draft" reads_1.fastq.gz reads_2.fastq.gz | samtools view --threads 3 -F 0x4 -b - | samtools fixmate -m --threads 3  - - | samtools sort -m 2g --threads 5 - | samtools markdup --threads 5 -r - sgs.sort.bam
        samtools index sgs.sort.bam
        python "$nextpolish_dir"/lib/nextpolish1.py -g "$draft" -t 1 -p 24 -s sgs.sort.bam > nextpolish_temp.fasta 2> nextpolish.txt
        bwa index nextpolish_temp.fasta
        bwa mem -t 24 nextpolish_temp.fasta reads_1.fastq.gz reads_2.fastq.gz | samtools view --threads 3 -F 0x4 -b - | samtools fixmate -m --threads 3  - - | samtools sort -m 2g --threads 5 - | samtools markdup --threads 5 -r - sgs.sort.bam
        samtools index sgs.sort.bam
        samtools faidx nextpolish_temp.fasta
        python "$nextpolish_dir"/lib/nextpolish1.py -g nextpolish_temp.fasta -t 2 -p 24 -s sgs.sort.bam 2>> nextpolish.txt | seqkit sort -l -r | seqtk seq -U > nextpolish.fasta
        rm sgs.sort.bam* nextpolish_temp.fasta*

        gzip nextpolish.fasta
    done
done
```

Pilon v1.24:
```bash
for d in $(seq -f "%04.1f" 0.1 0.1 50); do 
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        draft="$base_dir"/drafts_low_quality/"$a".fasta

        bwa mem -t 24 "$draft" reads_1.fastq.gz reads_2.fastq.gz | samtools sort > pilon_alignments.bam; samtools index pilon_alignments.bam
        pilon --genome "$draft" --frags pilon_alignments.bam --output pilon > pilon.txt 2>&1
        seqtk seq -U pilon.fasta > pilon_temp.fasta && mv pilon_temp.fasta pilon.fasta
        rm pilon_alignments.bam pilon_alignments.bam.bai

        gzip pilon.fasta
    done
done
```




## Process results

Assess assemblies with [`compare_assemblies.py`](https://github.com/rrwick/Perfect-bacterial-genome-tutorial/wiki/Comparing-assemblies) script:
```bash
for d in $(seq -f "%04.1f" 0.1 0.1 50); do 
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        ref="$base_dir"/references/"$a".fasta
        if [[ -f "polypolish.fasta.gz" && ! -f "polypolish.errors" ]]; then compare_assemblies.py --aligner edlib "$ref" polypolish.fasta.gz > polypolish.errors &; fi
        if [[ -f "polypolish-careful.fasta.gz" && ! -f "polypolish-careful.errors" ]]; then compare_assemblies.py --aligner edlib "$ref" polypolish-careful.fasta.gz > polypolish-careful.errors &; fi
        if [[ -f "pypolca.fasta.gz" && ! -f "pypolca.errors" ]]; then compare_assemblies.py --aligner edlib "$ref" pypolca.fasta.gz > pypolca.errors &; fi
        if [[ -f "pypolca-careful.fasta.gz" && ! -f "pypolca-careful.errors" ]]; then compare_assemblies.py --aligner edlib "$ref" pypolca-careful.fasta.gz > pypolca-careful.errors &; fi
        if [[ -f "hypo.fasta.gz" && ! -f "hypo.errors" ]]; then compare_assemblies.py --aligner edlib "$ref" hypo.fasta.gz > hypo.errors &; fi
        if [[ -f "fmlrc2.fasta.gz" && ! -f "fmlrc2.errors" ]]; then compare_assemblies.py --aligner edlib "$ref" fmlrc2.fasta.gz > fmlrc2.errors &; fi
        if [[ -f "nextpolish.fasta.gz" && ! -f "nextpolish.errors" ]]; then compare_assemblies.py --aligner edlib "$ref" nextpolish.fasta.gz > nextpolish.errors &; fi
        if [[ -f "pilon.fasta.gz" && ! -f "pilon.errors" ]]; then compare_assemblies.py --aligner edlib "$ref" pilon.fasta.gz > pilon.errors &; fi
    done
    wait
done
```

Produce TSV file of results:
```bash
printf "genome\ttarget_illumina_depth\tpolishing\tremaining_errors\n" > "$base_dir"/results_low_quality.tsv
for d in $(seq -f "%04.1f" 0.1 0.1 50); do 
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        ref="$base_dir"/references/"$a".fasta
        if [[ -f "polypolish.errors" ]]; then printf "$a\t$d\tPolypolish\t"$(cat polypolish.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results_low_quality.tsv; fi
        if [[ -f "polypolish-careful.errors" ]]; then printf "$a\t$d\tPolypolish-careful\t"$(cat polypolish-careful.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results_low_quality.tsv; fi
        if [[ -f "pypolca.errors" ]]; then printf "$a\t$d\tpypolca\t"$(cat pypolca.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results_low_quality.tsv; fi
        if [[ -f "pypolca-careful.errors" ]]; then printf "$a\t$d\tpypolca-careful\t"$(cat pypolca-careful.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results_low_quality.tsv; fi
        if [[ -f "hypo.errors" ]]; then printf "$a\t$d\tHyPo\t"$(cat hypo.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results_low_quality.tsv; fi
        if [[ -f "fmlrc2.errors" ]]; then printf "$a\t$d\tFMLRC2\t"$(cat fmlrc2.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results_low_quality.tsv; fi
        if [[ -f "nextpolish.errors" ]]; then printf "$a\t$d\tNextPolish\t"$(cat nextpolish.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results_low_quality.tsv; fi
        if [[ -f "pilon.errors" ]]; then printf "$a\t$d\tPilon\t"$(cat pilon.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results_low_quality.tsv; fi
    done
done
```



## Conclusions

* When draft genome error rates are very high, none of the tools were able to fix all of the errors. The Polypolish paper found that a multi-tool approach is best in a situation like this.
* FMLRC2 did the best at higher depths, consistent with the results from the FMLRC2 paper (10.1093/molbev/msad048). This is probably because it's not alignment based, which makes it better at fixing errors in repeats (see Figure S1).
* Compared to the main analysis using high quality draft genomes, the effect of introduced errors is much less obvious here. E.g. a few introduced errors isn't apparent when there are thousands of total errors.
* FMLRC2 and NextPolish still made the total error count worse at very low depths.
* Polypolish-careful did worse than Polypolish-defaults at higher depths because it struggles to fix errors in repeats.
