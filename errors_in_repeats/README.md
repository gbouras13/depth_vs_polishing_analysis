This test puts errors in a repeat sequence in the _Salmonella enterica_ genome and then runs each polishing tool using the 50x reads to see how well they can polish errors in repeats.



## Prepare test genome

The `ATCC_10708_Salmonella_enterica` genome has seven copies of the rRNA operon at these positions:
* 463686-468837
* 1130116-1135167
* 3652823-3658228
* 4402737-4407788
* 4445493-4450703
* 4600890-4605941
* 4695820-4701030


This Python code takes the reference genome and adds 5 substitution errors in each instance of the rRNA operon:
```python
import gzip
import pathlib
import random

random.seed(0)

operon_positions = [(463685, 468837), (1130115, 1135167), (3652822, 3658239),
                    (4402736, 4407788), (4445492, 4450703), (4600889, 4605941),
                    (4695819, 4701030)]
errors_per_operon = 5


def get_random_base():
    return {0: 'A', 1: 'C', 2: 'G', 3: 'T'}[random.randint(0, 3)]


def get_random_different_base(b):
    random_base = get_random_base()
    while b == random_base:
        random_base = get_random_base()
    return random_base


def add_errors(seq, error_count, start):
    seq = list(seq)
    for i in random.sample(range(len(seq)), error_count):
        print(i + start + 1)
        seq[i] = get_random_different_base(seq[i])
    return ''.join(seq)


reference_genome = pathlib.Path.cwd().parents[0] / 'reference_assemblies' / 'ATCC_10708_Salmonella_enterica.fasta.gz'

with gzip.open(reference_genome, 'rt') as f:
    lines = [line.rstrip() for line in f]

seq = lines[1]

for start, end in operon_positions:
    print(seq[start:end])
    # stay away from the edges of the repeat
    start += 1000
    end -= 1000
    repeat_seq = seq[start:end]
    repeat_seq_with_errors = add_errors(repeat_seq, errors_per_operon, start)
    seq = seq[:start] + repeat_seq_with_errors + seq[end:]

lines[1] = seq

with open('draft.fasta', 'wt') as f:
    for line in lines:
        f.write(line)
        f.write('\n')
```

The errors were introduced at these positions (1-based coordinates):
* 464851, 465746, 466263, 466408, 467790
* 1131504, 1131686, 1131688, 1132270, 1133183
* 3654235, 3655272, 3655601, 3655756, 3656115
* 4404736, 4405101, 4405758, 4406243, 4406728
* 4446822, 4446866, 4447803, 4448327, 4448717
* 4602265, 4603188, 4603466, 4603712, 4604332
* 4697102, 4697187, 4697353, 4699600, 4699923



## Run polishers

Set some variables:
```bash
base_dir=/home/wickr/2024-01_low-depth_Illumina_polishing
nextpolish_dir=/home/wickr/programs/NextPolish
genome_size=4801704
draft=draft.fasta
ref="$base_dir"/references/ATCC_10708_Salmonella_enterica.fasta
integer_depth=50
```

Prepare the directory:
```bash
mkdir "$base_dir"/errors_in_repeats
cd "$base_dir"/errors_in_repeats
cp "$base_dir"/50.0/ATCC_10708_Salmonella_enterica/reads_[12].fastq.gz .

# I then put the error-containing genome (created by the above Python code) into this directory

bwa index "$draft"
samtools faidx "$draft"
```

Polypolish v0.6.0 (both defaults and --careful):
```bash
bwa mem -t 24 -a "$draft" reads_1.fastq.gz > alignments_1.sam 2>> bwa.txt
bwa mem -t 24 -a "$draft" reads_2.fastq.gz > alignments_2.sam 2>> bwa.txt
polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam 2> polypolish_filter.txt
polypolish polish "$draft" filtered_1.sam filtered_2.sam 1> polypolish.fasta 2> polypolish.txt
polypolish polish --careful "$draft" filtered_1.sam filtered_2.sam 1> polypolish-careful.fasta 2> polypolish-careful.txt
rm *.sam
```

pypolca v0.3.0 (both defaults and --careful):
```bash
pypolca run -a "$draft" -1 reads_1.fastq.gz -2 reads_2.fastq.gz -t 24 -o pypolca > pypolca.txt 2>&1
seqtk seq pypolca/polca_corrected.fasta > pypolca.fasta
rm -r pypolca

pypolca run --careful -a "$draft" -1 reads_1.fastq.gz -2 reads_2.fastq.gz -t 24 -o pypolca > pypolca-careful.txt 2>&1
seqtk seq pypolca/polca_corrected.fasta > pypolca-careful.fasta
rm -r pypolca
```

HyPo v1.0.3:
```bash
bwa mem -t 24 "$draft" reads_1.fastq.gz reads_2.fastq.gz 2> /dev/null | samtools sort > alignments.bam; samtools index alignments.bam
echo -e "reads_1.fastq.gz\nreads_2.fastq.gz" > read_filenames.txt
hypo -d "$draft" -r @read_filenames.txt -s "$genome_size" -c "$integer_depth" -b alignments.bam -t 24 -o hypo.fasta > hypo.txt 2>&1
rm alignments.bam alignments.bam.bai read_filenames.txt
rm -r aux
```

FMLRC2 v0.1.8:
```bash
gunzip -c reads_1.fastq.gz reads_2.fastq.gz | awk 'NR % 4 == 2' | tr NT TN | ropebwt2 -LR | tr NT TN | fmlrc2-convert comp_msbwt.npy 2> fmlrc2.txt
fmlrc2 -t 24 comp_msbwt.npy "$draft" fmlrc2.fasta 2>> fmlrc2.txt
rm comp_msbwt.npy
```

NextPolish v1.4.1 (one round of their [user-defined short-read pipeline](https://nextpolish.readthedocs.io/en/latest/TUTORIAL.html)):
```bash
bwa mem -t 24 "$draft" reads_1.fastq.gz reads_2.fastq.gz | samtools view --threads 3 -F 0x4 -b - | samtools fixmate -m --threads 3  - - | samtools sort -m 2g --threads 5 - | samtools markdup --threads 5 -r - sgs.sort.bam
samtools index sgs.sort.bam
python "$nextpolish_dir"/lib/nextpolish1.py -g "$draft" -t 1 -p 24 -s sgs.sort.bam > nextpolish_temp.fasta 2> nextpolish.txt
bwa index nextpolish_temp.fasta
bwa mem -t 24 nextpolish_temp.fasta reads_1.fastq.gz reads_2.fastq.gz | samtools view --threads 3 -F 0x4 -b - | samtools fixmate -m --threads 3  - - | samtools sort -m 2g --threads 5 - | samtools markdup --threads 5 -r - sgs.sort.bam
samtools index sgs.sort.bam
samtools faidx nextpolish_temp.fasta
python "$nextpolish_dir"/lib/nextpolish1.py -g nextpolish_temp.fasta -t 2 -p 24 -s sgs.sort.bam 2>> nextpolish.txt | seqkit sort -l -r | seqtk seq -U > nextpolish.fasta 
rm sgs.sort.bam* nextpolish_temp.fasta*
```

Pilon v1.24:
```bash
bwa mem -t 24 "$draft" reads_1.fastq.gz reads_2.fastq.gz | samtools sort > alignments.bam; samtools index alignments.bam
pilon --genome "$draft" --frags alignments.bam --output pilon
seqtk seq -U pilon.fasta > temp.fasta && mv temp.fasta pilon.fasta
rm alignments.bam alignments.bam.bai
```



## Assess results

```bash
for p in polypolish polypolish-careful pypolca pypolca-careful fmlrc2 hypo nextpolish pilon; do
    compare_assemblies.py --aligner edlib "$ref" "$p".fasta > "$p".errors
done
```

```bash
printf "polishing_method\tremaining_errors\n" > results.tsv
for p in polypolish polypolish-careful pypolca pypolca-careful fmlrc2 hypo nextpolish pilon; do
    printf "$p\t"$(cat "$p".errors | grep -o "*" | wc -l)"\n" >> results.tsv
done
```

Error counts:
```
polypolish           1
polypolish-careful  35
pypolca             33
pypolca-careful     33
fmlrc2               0
hypo                35
nextpolish          30
pilon               32
```
