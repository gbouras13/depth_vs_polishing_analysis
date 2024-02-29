# Create Figure 1

* Took the first 151 bases of the ATCC 25922 E coli plasmid 4 in `ATCC_25922_plasmid_4.fasta`

* Changed base 49 T->G in `ATCC_25922_plasmid_4_with_error.fasta`


## High Confidence Error - Where Pypolca will always polish

* Ran [InSilicoSeq](https://doi.org/10.1093/bioinformatics/bty630) version 1.5.4

```
mamba install insilicoseq
iss generate --genomes ATCC_25922_plasmid_4_with_error.fasta --cpus 8 --model HiSeq  --compress --n_reads 10 --output ATCC_25922_plasmid_4_high_conf
```
* Then map the simulated reads to the reference with an error, sort the BAM file

```
minimap2 -ax sr ATCC_25922_plasmid_4.fasta ATCC_25922_plasmid_4_high_conf_R1.fastq.gz ATCC_25922_plasmid_4_high_conf_R2.fastq.gz > high_conf_alignment.sam
samtools view -bS high_conf_alignment.sam > high_conf_alignment.bam
samtools sort high_conf_alignment.bam -o high_conf_alignment_sorted.bam
samtools index high_conf_alignment_sorted.bam
```

## Pypolca Careful vs Pypolca Default low counts - Fig B

* To generate a mixture

```
iss generate --genomes ATCC_25922_plasmid_4_with_error.fasta --cpus 8 --model HiSeq  --compress --n_reads 5 --output ATCC_25922_plasmid_4_error_diff
iss generate --genomes ATCC_25922_plasmid_4.fasta --cpus 8 --model HiSeq  --compress --n_reads 2 --output ATCC_25922_plasmid_4_correct_diff

cat ATCC_25922_plasmid_4_error_diff_R1.fastq.gz ATCC_25922_plasmid_4_correct_diff_R1.fastq.gz > ATCC_25922_plasmid_4_combined_diff_R1.fastq.gz
cat ATCC_25922_plasmid_4_error_diff_R2.fastq.gz ATCC_25922_plasmid_4_correct_diff_R2.fastq.gz > ATCC_25922_plasmid_4_combined_diff_R2.fastq.gz
```
* Then map the simulated reads to the reference with an error, sort the BAM file

```
minimap2 -ax sr ATCC_25922_plasmid_4_with_error.fasta ATCC_25922_plasmid_4_combined_diff_R1.fastq.gz ATCC_25922_plasmid_4_combined_diff_R2.fastq.gz > diff_alignment.sam
samtools view -bS diff_alignment.sam > diff_alignment.bam
samtools sort diff_alignment.bam -o diff_alignment_sorted.bam
samtools index diff_alignment_sorted.bam
```

## Pypolca Careful vs Pypolca Default Ratio 

* To generate a mixture where ratio is between 2 and 3

```
iss generate --genomes ATCC_25922_plasmid_4_with_error.fasta --cpus 8 --model HiSeq  --compress --n_reads 14 --output ATCC_25922_plasmid_4_error_ratio
iss generate --genomes ATCC_25922_plasmid_4.fasta --cpus 8 --model HiSeq  --compress --n_reads 5 --output ATCC_25922_plasmid_4_correct_ratio

cat ATCC_25922_plasmid_4_error_ratio_R1.fastq.gz ATCC_25922_plasmid_4_correct_ratio_R1.fastq.gz > ATCC_25922_plasmid_4_combined_ratio_R1.fastq.gz
cat ATCC_25922_plasmid_4_error_ratio_R2.fastq.gz ATCC_25922_plasmid_4_correct_ratio_R2.fastq.gz > ATCC_25922_plasmid_4_combined_ratio_R2.fastq.gz
```
* Then map the simulated reads to the reference with an error, sort the BAM file

```
minimap2 -ax sr ATCC_25922_plasmid_4_with_error.fasta ATCC_25922_plasmid_4_combined_ratio_R1.fastq.gz ATCC_25922_plasmid_4_combined_ratio_R2.fastq.gz > ratio_alignment.sam
samtools view -bS ratio_alignment.sam > ratio_alignment.bam
samtools sort ratio_alignment.bam -o ratio_alignment_sorted.bam
samtools index ratio_alignment_sorted.bam
```