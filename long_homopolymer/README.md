This directory contains some analysis of the toughest single error to fix: an extra-long homopolymer in the _Salmonella enterica_ plasmid. 

The ONT-only assembly contained Cx20, while the reference sequence has Cx21:
```
ONT-only:  ...ACTGGCCCCCCCCCCCCCCCCCCCCGGATA...
reference: ...ACTGGCCCCCCCCCCCCCCCCCCCCCGGATA...
```



# Read counts

This code counts the exact-match instances of each length of the homopolymer in the raw reads, requiring 5 bp of adjacent sequence.
```bash
cd ~/2024-01_low-depth_Illumina_polishing

read_dir=/home/wickr/2023-10_PulseNet_ATCC_assemblies
r1="$read_dir"/ATCC_10708_Salmonella_enterica/reads/illumina_1.fastq.gz
r2="$read_dir"/ATCC_10708_Salmonella_enterica/reads/illumina_2.fastq.gz
ont="$read_dir"/ATCC_10708_Salmonella_enterica/reads/nanopore.fastq.gz

zcat "$r1" | grep -B1 -A2 -P --no-group-separator "ACTGGCCCCC+GGATA" > homopolymer_reads_illumina.fastq
zcat "$r1" | grep -B1 -A2 -P --no-group-separator "TATCCGGGGG+CCAGT" >> homopolymer_reads_illumina.fastq
zcat "$r2" | grep -B1 -A2 -P --no-group-separator "ACTGGCCCCC+GGATA" >> homopolymer_reads_illumina.fastq
zcat "$r2" | grep -B1 -A2 -P --no-group-separator "TATCCGGGGG+CCAGT" >> homopolymer_reads_illumina.fastq
zcat "$ont" | grep -B1 -A2 -P --no-group-separator "ACTGGCCCCC+GGATA" > homopolymer_reads_ont.fastq
zcat "$ont" | grep -B1 -A2 -P --no-group-separator "TATCCGGGGG+CCAGT" >> homopolymer_reads_ont.fastq

printf "len\tplatform\tbase\tcount\n" > homopolymer_counts.tsv
for len in $(seq 10 1 30); do
    c_homopolymer=$(printf '%.sC' {1.."$len"})
    g_homopolymer=$(printf '%.sG' {1.."$len"})
    c_with_neighbours="ACTGG"$c_homopolymer"GGATA"
    g_with_neighbours="TATCC"$g_homopolymer"CCAGT"
    illumina_c=$(cat homopolymer_reads_illumina.fastq | grep -c "$c_with_neighbours")
    illumina_g=$(cat homopolymer_reads_illumina.fastq | grep -c "$g_with_neighbours")
    ont_c=$(cat homopolymer_reads_ont.fastq | grep -c "$c_with_neighbours")
    ont_g=$(cat homopolymer_reads_ont.fastq | grep -c "$g_with_neighbours")
    printf "$len\tillumina\tc\t$illumina_c\n" >> homopolymer_counts.tsv
    printf "$len\tillumina\tg\t$illumina_g\n" >> homopolymer_counts.tsv
    printf "$len\tont\tc\t$ont_c\n" >> homopolymer_counts.tsv
    printf "$len\tont\tg\t$ont_g\n" >> homopolymer_counts.tsv
done
```



# Polisher performance


Check which polishers do best on that homopolymer:
```bash
zgrep -o "ACTGGCCCCCCCCCCCCCCCCCCCCCGGATA" */ATCC_10708_Salmonella_enterica/*.fasta.gz
```

Out of the 500 total assemblies for each polisher, this is how many had the correct homopolymer:
```
pypolca-defaults:     35
pypolca-careful:       2
polypolish-defaults:   0
polypolish-careful:    0
FMLRC2:                0
HyPo:                 19
NextPolish:          159
Pilon:                 3
```



# IGV

Prepare alignments for IGV using the full Illumina read set:
```bash
base_dir=/home/wickr/2024-01_low-depth_Illumina_polishing
read_dir=/home/wickr/2023-10_PulseNet_ATCC_assemblies
draft="$base_dir"/drafts/ATCC_10708_Salmonella_enterica.fasta
r1="$read_dir"/ATCC_10708_Salmonella_enterica/reads_qc/illumina_1.fastq.gz
r2="$read_dir"/ATCC_10708_Salmonella_enterica/reads_qc/illumina_2.fastq.gz

bwa mem -SP -t 24 "$draft" "$r1" "$r2" | samtools sort > alignments.bam
samtools index alignments.bam
```
