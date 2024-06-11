## Starting files

For each genome, I'm starting with the draft (ONT-only) and reference (manually curated ground truth):
```
drafts/ATCC_10708_Salmonella_enterica.fasta
drafts/ATCC_14035_Vibrio_cholerae.fasta
drafts/ATCC_17802_Vibrio_parahaemolyticus.fasta
drafts/ATCC_19119_Listeria_ivanovii.fasta
drafts/ATCC_25922_Escherichia_coli.fasta
drafts/ATCC_33560_Campylobacter_jejuni.fasta
drafts/ATCC_35221_Campylobacter_lari.fasta
drafts/ATCC_35897_Listeria_welshimeri.fasta
drafts/ATCC_BAA-679_Listeria_monocytogenes.fasta

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



## Illumina read QC

For each Illumina read set, I ran fastp like this:
```bash
for a in "${genomes[@]}"; do
    cd "$read_dir"/"$a"
    mkdir reads_qc
    fastp --in1 reads/illumina_1.fastq.gz --in2 reads/illumina_2.fastq.gz --out1 reads_qc/illumina_1.fastq.gz --out2 reads_qc/illumina_2.fastq.gz --unpaired1 reads_qc/illumina_u.fastq.gz --unpaired2 reads_qc/illumina_u.fastq.gz
done
```

I then used the `reads_qc/illumina_1.fastq.gz` and `reads_qc/illumina_2.fastq.gz` files for polishing and discarded the small `reads_qc/illumina_u.fastq.gz` file.



## Check Illumina depth for each genome

For each genome, count the total number of Illumina bases (post-QC) and divide by the genome size to get the depth:
```bash
for a in "${genomes[@]}"; do
    genome_size=${genome_sizes["$a"]}
    bases_1=$(fast_count "$read_dir"/"$a"/reads_qc/illumina_1.fastq.gz | cut -f3)
    bases_2=$(fast_count "$read_dir"/"$a"/reads_qc/illumina_2.fastq.gz | cut -f3)
    base_count=$(( $bases_1 + $bases_2 ))
    full_depth=$(echo "scale=6; $base_count / "$genome_size | bc)
    printf "$a\t$full_depth\n"
done
```
Each genome has >50x depth, so I'll test depths up to 50x.



## Build indices

Build BWA and samtools indices for the draft genomes:
```bash
cd "$base_dir"/drafts
for a in "${genomes[@]}"; do
    bwa index "$a".fasta
    samtools faidx "$a".fasta
done
```



## Subsample reads

I divided read bases by read count to get a mean read length for each genome:
```bash
declare -A mean_read_lengths
mean_read_lengths["ATCC_10708_Salmonella_enterica"]=112.4235209
mean_read_lengths["ATCC_14035_Vibrio_cholerae"]=127.3291922
mean_read_lengths["ATCC_17802_Vibrio_parahaemolyticus"]=121.8790019
mean_read_lengths["ATCC_19119_Listeria_ivanovii"]=105.0076326
mean_read_lengths["ATCC_25922_Escherichia_coli"]=121.9587802
mean_read_lengths["ATCC_33560_Campylobacter_jejuni"]=120.0284425
mean_read_lengths["ATCC_35221_Campylobacter_lari"]=128.8075044
mean_read_lengths["ATCC_35897_Listeria_welshimeri"]=129.3498371
mean_read_lengths["ATCC_BAA-679_Listeria_monocytogenes"]=120.3243788
```

Create the subsampled reads at 0.1x intervals up to 50.0x using [`seqtk sample`](https://github.com/lh3/seqtk). I used the read count as the random seed value so I could regenerate the exact same read sets as needed:
```bash
for d in $(seq -f "%04.1f" 0.1 0.1 50); do
    for a in "${genomes[@]}"; do
        mkdir -p "$base_dir"/"$d"/"$a"
        cd "$base_dir"/"$d"/"$a"
        echo "$a $d"x
        r1="$read_dir"/"$a"/reads_qc/illumina_1.fastq.gz
        r2="$read_dir"/"$a"/reads_qc/illumina_2.fastq.gz
        genome_size=${genome_sizes["$a"]}
        mean_read_length=${mean_read_lengths["$a"]}
        target_bases=$(echo "scale=2; $genome_size * "$d | bc)
        target_reads=$(echo "scale=2; $target_bases / "$mean_read_length" / 2" | bc | xargs printf "%1.0f")
        seqtk sample -s "$target_reads" "$r1" "$target_reads" > reads_1.fastq &
        seqtk sample -s "$target_reads" "$r2" "$target_reads" > reads_2.fastq &
        wait
        fast_count reads_1.fastq > reads_1.tsv
        fast_count reads_2.fastq > reads_2.tsv
    done
done
```

Gzip the reads to save space:
```bash
for d in $(seq -f "%04.1f" 0.1 0.1 50); do
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        pigz -p8 reads_1.fastq reads_2.fastq
    done
done
```

Make a read info table:
```bash
printf "genome\tgenome_size\ttarget_illumina_depth\tread_count\tread_bases\tactual_illumina_depth\n" > "$base_dir"/reads.tsv
for d in $(seq -f "%04.1f" 0.1 0.1 50); do
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        genome_size=${genome_sizes["$a"]}
        actual_read_count=$(( $(cut -f2 reads_1.tsv) + $(cut -f2 reads_2.tsv) ))
        actual_base_count=$(( $(cut -f3 reads_1.tsv) + $(cut -f3 reads_2.tsv) ))
        actual_depth=$(echo "scale=6; $actual_base_count / "$genome_size | bc)
        printf "$a\t$genome_size\t$d\t$actual_read_count\t$actual_base_count\t$actual_depth\n" >> "$base_dir"/reads.tsv
    done
done
```



## Run individual polishers

Polypolish v0.6.0 (both defaults and --careful):
```bash
for d in $(seq -f "%04.1f" 0.1 0.1 50); do
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        draft="$base_dir"/drafts/"$a".fasta

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
for d in $(seq -f "%04.1f" 0.1 0.1 50); do
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        draft="$base_dir"/drafts/"$a".fasta

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
        while [[ ! -f "hypo.fasta.gz" ]]; do  # HyPo randomly fails a lot when depth <5, so repeated tries are necessary.
            draft="$base_dir"/drafts/"$a".fasta
            genome_size=${genome_sizes["$a"]}
            integer_depth=$(printf "%.0f" "$d"); if [[ $integer_depth == 0 ]]; then integer_depth=1; fi

            bwa mem -t 24 "$draft" reads_1.fastq.gz reads_2.fastq.gz 2> /dev/null | samtools sort > alignments.bam; samtools index alignments.bam
            echo -e "reads_1.fastq.gz\nreads_2.fastq.gz" > read_filenames.txt
            hypo -d "$draft" -r @read_filenames.txt -s "$genome_size" -c "$integer_depth" -b alignments.bam -t 24 -o hypo.fasta > hypo.txt 2>&1
            rm alignments.bam alignments.bam.bai read_filenames.txt
            rm -r aux

            gzip hypo.fasta
        done
    done
done
```

FMLRC2 v0.1.8:
```bash
for d in $(seq -f "%04.1f" 0.1 0.1 50); do
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        draft="$base_dir"/drafts/"$a".fasta

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
        draft="$base_dir"/drafts/"$a".fasta

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
        draft="$base_dir"/drafts/"$a".fasta

        bwa mem -t 24 "$draft" reads_1.fastq.gz reads_2.fastq.gz | samtools sort > alignments.bam; samtools index alignments.bam
        pilon --genome "$draft" --frags alignments.bam --output pilon > pilon.txt 2>&1
        seqtk seq -U pilon.fasta > temp.fasta && mv temp.fasta pilon.fasta
        rm alignments.bam alignments.bam.bai

        gzip pilon.fasta
    done
done
```



## Run polisher combinations

Since Polypolish and pypolca were by far the best, I'm only trying combinations of those two.

All possible combinations of Polypolish and pypolca:
* Polypolish + pypolca
* Polypolish-careful + pypolca
* Polypolish + pypolca-careful
* Polypolish-careful + pypolca-careful
* pypolca + Polypolish
* pypolca-careful + Polypolish
* pypolca + Polypolish-careful
* pypolca-careful + Polypolish-careful

pypolca on Polypolish:
```bash
for d in $(seq -f "%04.1f" 0.1 0.1 50); do
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        draft="$base_dir"/drafts/"$a".fasta
        gunzip polypolish.fasta.gz polypolish-careful.fasta.gz

        pypolca run -a polypolish.fasta -1 reads_1.fastq.gz -2 reads_2.fastq.gz -t 24 -o pypolca > polypolish__pypolca.txt 2>&1
        seqtk seq pypolca/pypolca_corrected.fasta > polypolish__pypolca.fasta
        rm -r pypolca

        pypolca run --careful -a polypolish-careful.fasta -1 reads_1.fastq.gz -2 reads_2.fastq.gz -t 24 -o pypolca > polypolish-careful__pypolca-careful.txt 2>&1
        seqtk seq pypolca/pypolca_corrected.fasta > polypolish-careful__pypolca-careful.fasta
        rm -r pypolca

        pypolca run -a polypolish-careful.fasta -1 reads_1.fastq.gz -2 reads_2.fastq.gz -t 24 -o pypolca > polypolish-careful__pypolca.txt 2>&1
        seqtk seq pypolca/pypolca_corrected.fasta > polypolish-careful__pypolca.fasta
        rm -r pypolca

        pypolca run --careful -a polypolish.fasta -1 reads_1.fastq.gz -2 reads_2.fastq.gz -t 24 -o pypolca > polypolish__pypolca-careful.txt 2>&1
        seqtk seq pypolca/pypolca_corrected.fasta > polypolish__pypolca-careful.fasta
        rm -r pypolca

        gzip polypolish.fasta polypolish-careful.fasta polypolish__pypolca.fasta polypolish-careful__pypolca-careful.fasta polypolish-careful__pypolca.fasta polypolish__pypolca-careful.fasta
    done
done
```

Polypolish on pypolca:
```bash
for d in $(seq -f "%04.1f" 0.1 0.1 50); do
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        draft="$base_dir"/drafts/"$a".fasta
        gunzip pypolca.fasta.gz pypolca-careful.fasta.gz

        bwa index pypolca.fasta
        bwa mem -t 24 -a pypolca.fasta reads_1.fastq.gz > alignments_1.sam
        bwa mem -t 24 -a pypolca.fasta reads_2.fastq.gz > alignments_2.sam
        polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam 2> pypolca__polypolish_filter.txt
        polypolish polish pypolca.fasta filtered_1.sam filtered_2.sam 1> pypolca__polypolish.fasta 2> pypolca__polypolish.txt
        polypolish polish --careful pypolca.fasta filtered_1.sam filtered_2.sam 1> pypolca__polypolish-careful.fasta 2> pypolca__polypolish-careful.txt
        rm *.sam *.fasta.bwt *.fasta.pac *.fasta.ann *.fasta.amb *.fasta.sa

        bwa index pypolca-careful.fasta
        bwa mem -t 24 -a pypolca-careful.fasta reads_1.fastq.gz > alignments_1.sam
        bwa mem -t 24 -a pypolca-careful.fasta reads_2.fastq.gz > alignments_2.sam
        polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam 2> pypolca-careful__polypolish_filter.txt
        polypolish polish pypolca-careful.fasta filtered_1.sam filtered_2.sam 1> pypolca-careful__polypolish.fasta 2> pypolca-careful__polypolish.txt
        polypolish polish --careful pypolca-careful.fasta filtered_1.sam filtered_2.sam 1> pypolca-careful__polypolish-careful.fasta 2> pypolca-careful__polypolish-careful.txt
        rm *.sam *.fasta.bwt *.fasta.pac *.fasta.ann *.fasta.amb *.fasta.sa

        gzip pypolca.fasta pypolca-careful.fasta pypolca__polypolish.fasta pypolca-careful__polypolish-careful.fasta pypolca__polypolish-careful.fasta pypolca-careful__polypolish.fasta
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
        if [[ -f "polypolish__pypolca.fasta.gz" && ! -f "polypolish__pypolca.errors" ]]; then compare_assemblies.py --aligner edlib "$ref" polypolish__pypolca.fasta.gz > polypolish__pypolca.errors &; fi
        if [[ -f "polypolish-careful__pypolca-careful.fasta.gz" && ! -f "polypolish-careful__pypolca-careful.errors" ]]; then compare_assemblies.py --aligner edlib "$ref" polypolish-careful__pypolca-careful.fasta.gz > polypolish-careful__pypolca-careful.errors &; fi
        if [[ -f "polypolish-careful__pypolca.fasta.gz" && ! -f "polypolish-careful__pypolca.errors" ]]; then compare_assemblies.py --aligner edlib "$ref" polypolish-careful__pypolca.fasta.gz > polypolish-careful__pypolca.errors &; fi
        if [[ -f "polypolish__pypolca-careful.fasta.gz" && ! -f "polypolish__pypolca-careful.errors" ]]; then compare_assemblies.py --aligner edlib "$ref" polypolish__pypolca-careful.fasta.gz > polypolish__pypolca-careful.errors &; fi
        if [[ -f "pypolca__polypolish.fasta.gz" && ! -f "pypolca__polypolish.errors" ]]; then compare_assemblies.py --aligner edlib "$ref" pypolca__polypolish.fasta.gz > pypolca__polypolish.errors &; fi
        if [[ -f "pypolca-careful__polypolish-careful.fasta.gz" && ! -f "pypolca-careful__polypolish-careful.errors" ]]; then compare_assemblies.py --aligner edlib "$ref" pypolca-careful__polypolish-careful.fasta.gz > pypolca-careful__polypolish-careful.errors &; fi
        if [[ -f "pypolca__polypolish-careful.fasta.gz" && ! -f "pypolca__polypolish-careful.errors" ]]; then compare_assemblies.py --aligner edlib "$ref" pypolca__polypolish-careful.fasta.gz > pypolca__polypolish-careful.errors &; fi
        if [[ -f "pypolca-careful__polypolish.fasta.gz" && ! -f "pypolca-careful__polypolish.errors" ]]; then compare_assemblies.py --aligner edlib "$ref" pypolca-careful__polypolish.fasta.gz > pypolca-careful__polypolish.errors &; fi
    done
    wait
done
```

Produce TSV file of results:
```bash
printf "genome\ttarget_illumina_depth\tpolishing\tremaining_errors\n" > "$base_dir"/results.tsv
for d in $(seq -f "%04.1f" 0.1 0.1 50); do
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        ref="$base_dir"/references/"$a".fasta
        if [[ -f "polypolish.errors" ]]; then printf "$a\t$d\tPolypolish\t"$(cat polypolish.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
        if [[ -f "polypolish-careful.errors" ]]; then printf "$a\t$d\tPolypolish-careful\t"$(cat polypolish-careful.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
        if [[ -f "pypolca.errors" ]]; then printf "$a\t$d\tpypolca\t"$(cat pypolca.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
        if [[ -f "pypolca-careful.errors" ]]; then printf "$a\t$d\tpypolca-careful\t"$(cat pypolca-careful.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
        if [[ -f "hypo.errors" ]]; then printf "$a\t$d\tHyPo\t"$(cat hypo.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
        if [[ -f "fmlrc2.errors" ]]; then printf "$a\t$d\tFMLRC2\t"$(cat fmlrc2.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
        if [[ -f "nextpolish.errors" ]]; then printf "$a\t$d\tNextPolish\t"$(cat nextpolish.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
        if [[ -f "pilon.errors" ]]; then printf "$a\t$d\tPilon\t"$(cat pilon.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
        if [[ -f "polypolish__pypolca.errors" ]]; then printf "$a\t$d\tPolypolish+pypolca\t"$(cat polypolish__pypolca.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
        if [[ -f "polypolish-careful__pypolca-careful.errors" ]]; then printf "$a\t$d\tPolypolish-careful+pypolca-careful\t"$(cat polypolish-careful__pypolca-careful.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
        if [[ -f "polypolish-careful__pypolca.errors" ]]; then printf "$a\t$d\tPolypolish-careful+pypolca\t"$(cat polypolish-careful__pypolca.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
        if [[ -f "polypolish__pypolca-careful.errors" ]]; then printf "$a\t$d\tPolypolish+pypolca-careful\t"$(cat polypolish__pypolca-careful.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
        if [[ -f "pypolca__polypolish.errors" ]]; then printf "$a\t$d\tpypolca+Polypolish\t"$(cat pypolca__polypolish.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
        if [[ -f "pypolca-careful__polypolish-careful.errors" ]]; then printf "$a\t$d\tpypolca-careful+Polypolish-careful\t"$(cat pypolca-careful__polypolish-careful.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
        if [[ -f "pypolca__polypolish-careful.errors" ]]; then printf "$a\t$d\tpypolca+Polypolish-careful\t"$(cat pypolca__polypolish-careful.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
        if [[ -f "pypolca-careful__polypolish.errors" ]]; then printf "$a\t$d\tpypolca-careful+Polypolish\t"$(cat pypolca-careful__polypolish.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
    done
done
```



## Hapo-G

The paper didn't include [Hapo-G](https://github.com/institut-de-genomique/HAPO-G), so this analysis occurred after the paper's publication.

Run Hapo-G v1.3.8:
```bash
conda activate hapog

for d in $(seq -f "%04.1f" 0.1 0.1 50); do
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        draft="$base_dir"/drafts/"$a".fasta

        hapog --genome "$draft" --pe1 reads_1.fastq.gz --pe2 reads_2.fastq.gz -o hapog -t 24 -u > hapog.txt 2>&1
        seqtk seq hapog/hapog_results/hapog.fasta > hapog.fasta
        rm -r hapog

        gzip hapog.fasta
    done
done
```

Assess assemblies:
```bash
for d in $(seq -f "%04.1f" 0.1 0.1 50); do
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        ref="$base_dir"/references/"$a".fasta
        if [[ -f "hapog.fasta.gz" && ! -f "hapog.errors" ]]; then compare_assemblies.py --aligner edlib "$ref" hapog.fasta.gz > hapog.errors &; fi
    done
    wait
done
```

Add to results TSV:
```bash
for d in $(seq -f "%04.1f" 0.1 0.1 50); do
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        ref="$base_dir"/references/"$a".fasta
        if [[ -f "hapog.errors" ]]; then printf "$a\t$d\tHapo-G\t"$(cat hapog.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
    done
done
```
