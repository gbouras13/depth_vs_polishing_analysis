## Starting files

* For the `hybracter` analysis, we decided to run all assemblies from scratch implementing our recommended best practices in the `main_analysis` benchmarking. These changes were implemented in `hybracter` v0.7.0.
* This also takes into account various factors such as assembler non-determinism (most assemblers, including Flye, are [not deterministic](https://plassembler.readthedocs.io/en/latest/flye_non_determinism/) ).
* To test out the effect of short read polishing at different depths, we pre subsampled the short reads similar to the `main_analysis`
* We used `hybracter v0.7.0` with default settings (which will subsampling the long reads to the best estimated 100x coverage read set - in practice slightly less here depending on the `-c` minimum chromosome length used)
* We also conducted no QC prior to `hybracter`, as it is included by default.

* Note for the `ATCC_10708_Salmonella_enterica`, the long read FASTQ file was enormous (15GB) so I subsampled it first

```
mv SRR27638402.fastq.gz SRR27638402_original.fastq.gz
filtlong --target_bases 1000000000 --min_mean_q 15 --min_length 1000 SRR27638402_original.fastq.gz | pigz > SRR27638402.fastq.gz
```





## Subsample reads

Taking the mean read length for each genome analysis conducted in the `main_analysis`:

```bash
base_dir=/home/user/Documents/hybracter_polishing_analysis/short_reads_subsampled
read_dir=/home/user/Documents/hybracter_polishing_analysis/short_reads
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
for d in $(seq -f "%02.0f" 1 1 50); do 
    for a in "${genomes[@]}"; do
        mkdir -p "$base_dir"/"$d"/"$a"
        cd "$base_dir"/"$d"/"$a"
        echo "$a $d"x
        r1="$read_dir"/"$a"/illumina_1.fastq.gz
        r2="$read_dir"/"$a"/illumina_2.fastq.gz
        genome_size=${genome_sizes["$a"]}
        mean_read_length=${mean_read_lengths["$a"]}
        target_bases=$(echo "scale=2; $genome_size * "$d | bc)
        target_reads=$(echo "scale=2; $target_bases / "$mean_read_length" / 2" | bc | xargs printf "%1.0f")
        seqtk sample -s "$target_reads" "$r1" "$target_reads" > reads_1.fastq &
        seqtk sample -s "$target_reads" "$r2" "$target_reads" > reads_2.fastq &
        wait
        seqtk size reads_1.fastq > reads_1.tsv
        seqtk size reads_2.fastq > reads_2.tsv
    done
done

```

Gzip the reads to save space:
```bash
for d in $(seq -f "%02.0f" 1 1 50); do 
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        pigz -p8 reads_1.fastq reads_2.fastq
    done
done
```



Make a read info table:
```bash
printf "genome\tgenome_size\ttarget_illumina_depth\tread_count\tread_bases\tactual_illumina_depth\n" > "$base_dir"/reads.tsv
for d in $(seq -f "%02.0f" 1 1 50); do 
    for a in "${genomes[@]}"; do
        cd "$base_dir"/"$d"/"$a"
        genome_size=${genome_sizes["$a"]}
        actual_read_count=$(( $(cut -f1 reads_1.tsv) + $(cut -f1 reads_2.tsv) ))
        actual_base_count=$(( $(cut -f2 reads_1.tsv) + $(cut -f2 reads_2.tsv) ))
        actual_depth=$(echo "scale=6; $actual_base_count / "$genome_size | bc)
        printf "$a\t$genome_size\t$d\t$actual_read_count\t$actual_base_count\t$actual_depth\n" >> "$base_dir"/reads.tsv
    done
done
```

### hybracter

* per sample to save some space
* One QC step before to save time 
* Buy default, this was with `--logic last`


```
mkdir -p hybracter_outputs


hybracter hybrid -i csvs/ATCC_33560_Campylobacter_jejuni.csv -o hybracter_outputs/ATCC_33560_Campylobacter_jejuni --threads 32 --no_medaka --min_quality 15 --subsample_depth 100 --configfile bhinion_config.yaml 

hybracter hybrid -i csvs/ATCC_35221_Campylobacter_lari.csv -o hybracter_outputs/ATCC_35221_Campylobacter_lari --threads 32 --no_medaka --min_quality 15 --subsample_depth 100 --configfile bhinion_config.yaml  

hybracter hybrid -i csvs/ATCC_25922_Escherichia_coli.csv -o hybracter_outputs/ATCC_25922_Escherichia_coli --threads 32 --no_medaka --min_quality 15 --subsample_depth 100 --configfile bhinion_config.yaml  

hybracter hybrid -i csvs/ATCC_19119_Listeria_ivanovii.csv -o hybracter_outputs/ATCC_19119_Listeria_ivanovii --threads 32 --no_medaka --min_quality 15 --subsample_depth 100 --configfile bhinion_config.yaml  

hybracter hybrid -i csvs/ATCC_BAA-679_Listeria_monocytogenes.csv -o hybracter_outputs/ATCC_BAA-679_Listeria_monocytogenes --threads 32 --no_medaka --min_quality 15 --subsample_depth 100 --configfile bhinion_config.yaml  

hybracter hybrid -i csvs/ATCC_10708_Salmonella_enterica.csv -o hybracter_outputs/ATCC_10708_Salmonella_enterica --threads 32 --no_medaka --min_quality 15 --subsample_depth 100 --configfile bhinion_config.yaml  

hybracter hybrid -i csvs/ATCC_35897_Listeria_welshimeri.csv -o hybracter_outputs/ATCC_35897_Listeria_welshimeri --threads 32 --no_medaka --min_quality 15 --subsample_depth 100 --configfile bhinion_config.yaml  

# due to double chrom 1M at 400 is approx 100x
hybracter hybrid -i csvs/ATCC_14035_Vibrio_cholerae.csv -o hybracter_outputs/ATCC_14035_Vibrio_cholerae --threads 32 --no_medaka --min_quality 15 --subsample_depth 400 --configfile bhinion_config.yaml  

# due to double chrom 1.65M at 300 is approx 100x overall
hybracter hybrid -i csvs/ATCC_17802_Vibrio_parahaemolyticus.csv -o hybracter_outputs/ATCC_17802_Vibrio_parahaemolyticus --threads 32 --no_medaka --min_quality 15 --subsample_depth 300 --configfile bhinion_config.yaml  

```




```bash 

mamba create -n compare_assemblies mappy edlib
conda activate compare_assemblies
```




## Process results

* Assess assemblies with [`compare_assemblies.py`](https://github.com/rrwick/Perfect-bacterial-genome-tutorial/wiki/Comparing-assemblies) script:

* _Vibrio cholerae_ was excluded due to the known structural heterogeneity in assembly methods with hybracter.
* I focused only on chromosomes
* Therefore, I made reference chromosome assemblies in `reference_chromosome_assemblies_hybracter` - these removed the plasmids from _E coli_ and _S enterica_, as we wanted to compare the impact of polishing on the chromosome. This was due to the large number of false positive plasmid contigs generated (due to barcode bleed observed in the short read sets)
* Additionally, I reoriented _V. parahaemolyticus_ with [Dnaapler](https://github.com/gbouras13/dnaapler) v0.7.0 so that the smaller chromosome began in a consistent fashion with a homolog of the repA gene as in hybracter to enable comparison.

```
dnaapler all -i references_assemblies/ATCC_17802_Vibrio_parahaemolyticus.fasta -o ATCC_17802_Vibrio_parahaemolyticus_dnaapler -t 8
mv ATCC_17802_Vibrio_parahaemolyticus_dnaapler/dnaapler_reoriented.fasta    reference_chromosome_assemblies_hybracter/ATCC_17802_Vibrio_parahaemolyticus.fasta
```

* ALE often preferred worse pre-polished assemblies as best (aka `--logic best` in hybracter), so I decided to also check this with the `pypolca` output aka what you would get if you ran `-logic last`


```bash

mamba create -n compare_assemblies edlib mappy

base_dir=/home/user/Documents/hybracter_polishing_analysis
read_dir=/home/user/Documents/hybracter_polishing_analysis/short_reads
genomes=(
    "ATCC_10708_Salmonella_enterica"
#    "ATCC_14035_Vibrio_cholerae"
    "ATCC_17802_Vibrio_parahaemolyticus"
    "ATCC_19119_Listeria_ivanovii"
    "ATCC_25922_Escherichia_coli"
    "ATCC_33560_Campylobacter_jejuni"
    "ATCC_35221_Campylobacter_lari"
    "ATCC_35897_Listeria_welshimeri"
    "ATCC_BAA-679_Listeria_monocytogenes"
)

mkdir -p "$base_dir"/compare_assemblies
conda activate compare_assemblies
for d in {1..50}; do 
    for a in "${genomes[@]}"; do
        mkdir -p "$base_dir"/compare_assemblies/"$d"/"$a"
        cd "$base_dir"/compare_assemblies/"$d"/"$a"
        ref="$base_dir"/reference_chromosome_assemblies_hybracter/"$a".fasta
        best_assembly="$base_dir"/hybracter_outputs/"$a"/FINAL_OUTPUT/complete/"$d"x_chromosome.fasta
        "$base_dir"/compare_assemblies.py --aligner edlib "$ref" $best_assembly > hybracter_best.errors 
        last_assembly="$base_dir"/hybracter_outputs/"$a"/supplementary_results/intermediate_chromosome_assemblies/"$d"x/"$d"x_pypolca.fasta 
        "$base_dir"/compare_assemblies.py --aligner edlib "$ref" $last_assembly > hybracter_last.errors 
    done
    wait
done
```


Produce TSV file of results:
```bash

base_dir=/home/user/Documents/hybracter_polishing_analysis
read_dir=/home/user/Documents/hybracter_polishing_analysis/short_reads

printf "genome\ttarget_illumina_depth\tpolishing\tremaining_errors\n" > "$base_dir"/results.tsv
for d in {1..50}; do 
    for a in "${genomes[@]}"; do
        cd "$base_dir"/compare_assemblies/"$d"/"$a"
        ref="$base_dir"/reference_chromosome_assemblies_hybracter/"$a".fasta
        if [[ -f "hybracter_best.errors" ]]; then printf "$a\t$d\thybracter_best\t"$(cat hybracter_best.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
        if [[ -f "hybracter_last.errors" ]]; then printf "$a\t$d\thybracter_last\t"$(cat hybracter_last.errors | grep -o "*" | wc -l)"\n" >> "$base_dir"/results.tsv; fi
    done
done
```

Summary TSVs

```bash
mkdir -p hybracter_summary_tsvs

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


for a in "${genomes[@]}"; do
    base_dir=/home/user/Documents/hybracter_polishing_analysis

    cp "$base_dir"/hybracter_outputs/"$a"/FINAL_OUTPUT/hybracter_summary.tsv hybracter_summary_tsvs/"$a".tsv

done

```




