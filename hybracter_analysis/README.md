## Starting files

For each genome, I'm starting with the draft (ONT-only) and reference (manually curated ground truth):
```
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

* For the `hybracter` analysis, we decided to run all assemblies from scratch implementing our recommended best practices in the `main_analysis` benchmarking. These changes were implemented in `hybracter` v0.7.0.
* This also takes into account various factors such as assembler non-determinism (most assemblers, including Flye, are [not deterministic](https://plassembler.readthedocs.io/en/latest/flye_non_determinism/) ).
* To test out the effect of short read polishing at different depths, we pre subsampled the short reads similar to the `main_analysis`
* We used `hybracter v0.7.0` with default settings (which will subsampling the long reads to estimated 100x coverage)

## Subsample reads

Taking the mean read length for each genome analysis conducted in teh `main_analysis`:

```bash
read_dir=/home/wickr/2023-10_PulseNet_ATCC_assemblies
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

* Instead of subsampling at 0.1x intervals, we subsampled at 1x intervals starting from 1x, due to the fact `hybracter` will create a long read assembly for every isolate are therefore require a lot more resources than the main analysis.  The same code using [`seqtk sample`](https://github.com/lh3/seqtk) was used. I used the read count as the random seed value so I could regenerate the exact same read sets as needed:

```bash
for d in $(seq -f "%04.1f" 1 1 50); do 
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
