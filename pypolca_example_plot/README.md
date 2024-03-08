# Steps to Create Figure 1

* I decided the easiest way to do this was to simulate fake reads and repurpose some code from [other work](https://github.com/gbouras13/CRS_Saureus_Evolutionary_Landscape).

* I took the first 151 bases of the ATCC 25922 E coli plasmid 4 in `ATCC_25922_plasmid_4.fasta`
* I then manually changed base 49 T->G in `ATCC_25922_plasmid_4_with_error.fasta`


## Figure 1A - High Confidence Error - Where Pypolca will always polish

* Then I ran [InSilicoSeq](https://doi.org/10.1093/bioinformatics/bty630) version 1.5.4

```bash
mamba install insilicoseq
iss generate --genomes ATCC_25922_plasmid_4_with_error.fasta --cpus 8 --model HiSeq  --compress --n_reads 10 --output ATCC_25922_plasmid_4_high_conf
```
* The simulated reads were then mapped to the reference with an error, followed by sorting the BAM file

```bash
minimap2 -ax sr ATCC_25922_plasmid_4.fasta ATCC_25922_plasmid_4_high_conf_R1.fastq.gz ATCC_25922_plasmid_4_high_conf_R2.fastq.gz > high_conf_alignment.sam
samtools view -bS high_conf_alignment.sam > high_conf_alignment.bam
samtools sort high_conf_alignment.bam -o high_conf_alignment_sorted.bam
samtools index high_conf_alignment_sorted.bam
```

## Figure 1B - Pypolca Careful vs Pypolca Default low counts 

* To generate a mixture and combine

```bash
iss generate --genomes ATCC_25922_plasmid_4_with_error.fasta --cpus 8 --model HiSeq  --compress --n_reads 5 --output ATCC_25922_plasmid_4_error_diff
iss generate --genomes ATCC_25922_plasmid_4.fasta --cpus 8 --model HiSeq  --compress --n_reads 2 --output ATCC_25922_plasmid_4_correct_diff

cat ATCC_25922_plasmid_4_error_diff_R1.fastq.gz ATCC_25922_plasmid_4_correct_diff_R1.fastq.gz > ATCC_25922_plasmid_4_combined_diff_R1.fastq.gz
cat ATCC_25922_plasmid_4_error_diff_R2.fastq.gz ATCC_25922_plasmid_4_correct_diff_R2.fastq.gz > ATCC_25922_plasmid_4_combined_diff_R2.fastq.gz
```
* Then map the simulated reads to the reference with an error, sort the BAM file

```bash
minimap2 -ax sr ATCC_25922_plasmid_4_with_error.fasta ATCC_25922_plasmid_4_combined_diff_R1.fastq.gz ATCC_25922_plasmid_4_combined_diff_R2.fastq.gz > diff_alignment.sam
samtools view -bS diff_alignment.sam > diff_alignment.bam
samtools sort diff_alignment.bam -o diff_alignment_sorted.bam
samtools index diff_alignment_sorted.bam
```

## Figure 1C - Pypolca Careful vs Pypolca Default Ratio 

* To generate a mixture where ratio is between 2 and 3

```bash
iss generate --genomes ATCC_25922_plasmid_4_with_error.fasta --cpus 8 --model HiSeq  --compress --n_reads 14 --output ATCC_25922_plasmid_4_error_ratio
iss generate --genomes ATCC_25922_plasmid_4.fasta --cpus 8 --model HiSeq  --compress --n_reads 5 --output ATCC_25922_plasmid_4_correct_ratio

cat ATCC_25922_plasmid_4_error_ratio_R1.fastq.gz ATCC_25922_plasmid_4_correct_ratio_R1.fastq.gz > ATCC_25922_plasmid_4_combined_ratio_R1.fastq.gz
cat ATCC_25922_plasmid_4_error_ratio_R2.fastq.gz ATCC_25922_plasmid_4_correct_ratio_R2.fastq.gz > ATCC_25922_plasmid_4_combined_ratio_R2.fastq.gz
```
* Then map the simulated reads to the reference with an error, sort the BAM file

```bash
minimap2 -ax sr ATCC_25922_plasmid_4_with_error.fasta ATCC_25922_plasmid_4_combined_ratio_R1.fastq.gz ATCC_25922_plasmid_4_combined_ratio_R2.fastq.gz > ratio_alignment.sam
samtools view -bS ratio_alignment.sam > ratio_alignment.bam
samtools sort ratio_alignment.bam -o ratio_alignment_sorted.bam
samtools index ratio_alignment_sorted.bam
```

* Then run `create_gviz_plots.R` to make the 3 subfigures. I manually combined them for the paper.