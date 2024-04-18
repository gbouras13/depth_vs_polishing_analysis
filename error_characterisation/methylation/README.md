A deep-dive on methylation for ATCC 35221 _C. lari_ was done using [MicrobeMod](https://github.com/cultivarium/MicrobeMod)

* 4mC, 5mC and 6mA models were tested.

First, we ran Dorado with methylation aware basecalling

```
./dorado-0.6.0-linux-x64/bin/dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v4.3.0 C_lari_pod5s  --modified-bases-models dna_r10.4.1_e8.2_400bps_sup@v4.3.0_6mA@v2,res_dna_r10.4.1_e8.2_400bps_sup@v4.3.0_4mC_5mC@v1 > c_lari.bam
```

Then we ran MicrobeMod

```
samtools fastq c_lari.bam -T MM,ML | minimap2 -t 32 --secondary=no -ax map-ont -y ATCC_35221_Campylobacter_lari.fasta -| \
samtools view -b | samtools sort -@ 32 -o c_lari.mapped.bam

samtools index c_lari.mapped.bam

MicrobeMod call_methylation -b c_lari.mapped.bam -r ATCC_35221_Campylobacter_lari.fasta -t 32

MicrobeMod call_methylation -b c_lari.mapped.bam -r ATCC_35221_Campylobacter_lari.fasta -t 32 --percent_methylation_cutoff 0.5  --output_prefix c_lari_0.5
```

The relevant outputs are `c_lari.mapped_methylated_sites.tsv` which contained each methylated site, and `c_lari.mapped_motifs.tsv`, which provides information on the identified methylation motifs.
