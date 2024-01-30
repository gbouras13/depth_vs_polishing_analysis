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

zcat "$r1" | grep -P "ACTGGCCCCC+GGATA" > c_seqs_r1.txt
zcat "$r1" | grep -P "TATCCGGGGG+CCAGT" > g_seqs_r1.txt
zcat "$r2" | grep -P "ACTGGCCCCC+GGATA" > c_seqs_r2.txt
zcat "$r2" | grep -P "TATCCGGGGG+CCAGT" > g_seqs_r2.txt

printf "len\tc_count\tg_count\n" > homopolymer_counts.tsv
for len in $(seq 10 1 30); do
    c_homopolymer=$(printf '%.sC' {1.."$len"})
    g_homopolymer=$(printf '%.sG' {1.."$len"})
    c_with_neighbours="ACTGG"$c_homopolymer"GGATA"
    g_with_neighbours="TATCC"$g_homopolymer"CCAGT"
    c_count=$(cat c_seqs_r*.txt | grep -c "$c_with_neighbours")
    g_count=$(cat g_seqs_r*.txt | grep -c "$g_with_neighbours")
    printf "$len\t$c_count\t$g_count\n" >> homopolymer_counts.tsv
done

rm *_seqs_r*.txt
```

```r
homopolymer <- data.frame(len=  c(17, 18, 19, 20, 21, 22, 23),
                          count=c( 0,  3,  3,  9, 15,  2,  0))

ggplot(homopolymer, aes(x=len, y=count)) +
  geom_bar(stat="identity", fill="#1f78b4") +
  theme_bw() + theme(panel.grid.major.x = element_blank()) +
  scale_x_continuous(limits = c(17.5, 22.5), expand = c(0, 0), minor_breaks=NULL) +
  scale_y_continuous(limits = c(0, 16), expand = c(0, 0), minor_breaks=seq(0, 16, 1)) +
  labs(title=NULL, x="homopolymer length", y="read count")
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

NextPolish did best, followed by pypolca. None of the polishers were able to reliably get the correct answer.

```r
fixed_counts <- data.frame(polisher=c("pypolca-defaults", "pypolca-careful", "Polypolish-defaults", "Polypolish-careful", "FMLRC2", "HyPo", "NextPolish", "Pilon"),
                          correct=c(35, 2, 0, 0, 0, 19, 159, 3))
fixed_counts$incorrect <- 500 - fixed_counts$correct

fixed_counts %>% pivot_longer(!polisher, names_to = "result", values_to = "count") -> fixed_counts
fixed_counts$result <- factor(fixed_counts$result, levels = c("incorrect", "correct"))
fixed_counts$polisher <- factor(fixed_counts$polisher, levels = c("pypolca-defaults", "pypolca-careful", "Polypolish-defaults", "Polypolish-careful", "FMLRC2", "HyPo", "NextPolish", "Pilon"))

ggplot(fixed_counts, aes(x=polisher, y=count, fill=result)) +
  geom_bar(stat="identity") +
  theme_bw() + theme(legend.title=element_blank()) +
  scale_fill_manual(values = c("correct" = "#1f78b4",
                               "incorrect" = "#fc8d62")) +
  scale_y_continuous(limits = c(0, 500), expand = c(0, 0), minor_breaks=NULL) +
  scale_x_discrete(labels = c("pypolca-\ndefaults", "pypolca\n-careful", "Polypolish-\ndefaults", "Polypolish-\ncareful", "FMLRC2", "HyPo", "NextPolish", "Pilon")) +
  labs(title=NULL, x=NULL, y="assembly count")
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
