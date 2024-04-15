This analysis does a parameter sweep for each polisher, trying to find a set of parameters that does well with low-depth polishing. For this test, I used only the _Campylobacter lari_ genome, because it was the smallest genome (so fastest to run) and had the highest density of pre-polishing errors. I also just tested four depths (5x, 10x, 15x and 20x) to save time. Parameters were ranked on the total number of remaining errors across those four depths. Since the _Campylobacter lari_ genome had 18 errors, a polisher that did nothing would result in a total of 72 errors.

Each parameter sweep is done by a Python script. The script first tries the default parameters, and if the documentation/publication suggests any alternative parameters (e.g. 'eukaryote mode' in FMLRC2), then it tries those also. Then it alternates between trying a random set of parameters and a mutated set of parameters. Mutated parameters start with the best set since so far and adjust the values slightly. A total of 1000 parameter sets were tried for each tool.

For each tool, I did my best to decide which parameters were relevant and which weren't, and only search relevant parameters. For example, a minimum-depth parameter would be relevant but a thread-count parameter would not. I also did my best to test each parameter in a range which makes sense, e.g. a 'depth' parameter might range from 0-50, while a 'fraction' parameter might range from 0-1.

Often, a tool had multiple similar parameter sets which tied for the best. To get a single 'best' set, I did my best to take a median set. E.g. if the best sets had a parameter with values 0.21, 0.22 or 0.23, I would use the set with 0.22.



## Totals for defaults

Before doing the parameter sweep, I got these totals (error counts for _C. lari_ across the four tested depths) from the main analysis results.

* Polypolish-default:   31
* Polypolish-careful:   31
* Pypolca-default:      63
* Pypolca-careful:      17
* FMLRC2-default:     2514
* HyPo:                105
* NextPolish:          960
* Pilon:               778

This shows that using defaults, Polypolish and Pypolca were the only tools that improved the _C. lari_ assembly (totals of <72).



## Set variables

```bash
base_dir=/home/wickr/2024-01_low-depth_Illumina_polishing
read_dir=/home/wickr/2023-10_PulseNet_ATCC_assemblies
nextpolish_dir=/home/wickr/programs/NextPolish
```



## Create directories

```bash
base_dir=/home/wickr/2024-01_low-depth_Illumina_polishing
cd "$base_dir"
mkdir parameter_sweep

cd "$base_dir"/parameter_sweep
mkdir fmlrc2
mkdir polypolish
mkdir pypolca
mkdir hypo
mkdir nextpolish
mkdir pilon
```



## Choosing best parameters

Often many parameter sets will tie for the best, so this Python code lets me choose whichever one is closest to the defaults:
```python
import math
from statistics import mean, stdev

def normalise_parameters(parameters):
    means = [mean(position) for position in zip(*parameters)]
    std_devs = [stdev(position) for position in zip(*parameters)]
    normalized_parameters = [tuple((element - means[i]) / (std_devs[i] if std_devs[i] != 0 else 1) for i, element in enumerate(tup)) for tup in parameters]
    return normalized_parameters, means, std_devs

def euclidean_distance(tup1, tup2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(tup1, tup2)))

def find_closest_parameters(all_parameters, target_parameters):
    normalized_parameters, means, std_devs = normalise_parameters(all_parameters)
    normalized_target = tuple((element - means[i]) / (std_devs[i] if std_devs[i] != 0 else 1) for i, element in enumerate(target_parameters))
    distances = [euclidean_distance(normalized_target, tup) for tup in normalized_parameters]
    min_index = distances.index(min(distances))
    return all_parameters[min_index]
```



## FMLRC2

```bash
cd "$base_dir"/parameter_sweep/fmlrc2
conda activate assembly
./fmlrc2.py
rm comp_msbwt_* fmlrc2.fasta
cat results.tsv | sort -t$'\t' -nk8,8 | head -n2
```
Best parameters: `-k 27 -m 0 -f 0.29`
Best total: 26



## HyPo

```bash
cd "$base_dir"/parameter_sweep/hypo
conda activate hypo
./hypo.py
rm alignments_*.bam* read_filenames_*.txt hypo.fasta
cat results.tsv | sort -t$'\t' -nk8,8 | head -n2
```
Best parameters: `-m 5 -x -5 -g -8 -n 20 -q 41`
Best total: 102

Most important parameter was `-q` (read alignment quality threshold).



## NextPolish

```bash
cd "$base_dir"/parameter_sweep/nextpolish
conda activate nextpolish
./nextpolish.py
rm sgs.sort_*.bam*
cat results.tsv | sort -t$'\t' -nk8,8 | head -n2
```
Best parameters: `-min_map_quality 1 -max_ins_fold_sgs 7 -max_clip_ratio_sgs 0 -trim_len_edge 0 -ext_len_edge 0 -indel_balance_factor_sgs 0.63 -min_count_ratio_skip 0.55 -max_len_kmer 35 -min_len_inter_kmer 5 -max_count_kmer 47`
Best total: 596



## Pilon

```bash
cd "$base_dir"/parameter_sweep/pilon
conda activate pilon
./pilon.py
rm alignments_*.bam*
cat results.tsv | sort -t$'\t' -nk8,8 | head -n2
```
Best parameters: `--fix bases --flank 3 --K 51 --mindepth 3 --minmq 58 --minqual 15`
Best total: 19

Most important parameter was either turning off reassembly (`--fix bases`) or using a large k-mer for reassembly.



## Polypolish

```bash
cd "$base_dir"/parameter_sweep/polypolish
conda activate assembly
./polypolish.py
rm *.sam polypolish.fasta
cat results.tsv | sort -t$'\t' -nk8,8 | head -n2
```
Best parameters: `-i 0.42 -v 0.68 -m 3 -d 3`
Best total: 17



## Pypolca

```bash
cd "$base_dir"/parameter_sweep/pypolca
conda activate assembly
./pypolca.py
cat results.tsv | sort -t$'\t' -nk8,8 | head -n2
```
Best parameters: `--min_alt 3 --min_ratio 3.5`
Best total: 14



## Best total for each tool

* Polypolish:  17
* Pypolca:     14
* FMLRC2:      26
* HyPo:       102
* NextPolish: 596
* Pilon:       19



## Combine results into one table

```bash
cd "$base_dir"/parameter_sweep
cat fmlrc2/results.tsv > results.tsv
tail -n+2 hypo/results.tsv >> results.tsv
tail -n+2 nextpolish/results.tsv >> results.tsv
tail -n+2 pilon/results.tsv >> results.tsv
tail -n+2 polypolish/results.tsv >> results.tsv
tail -n+2 pypolca/results.tsv >> results.tsv
```



## Conclusions

* FMLRC2 and Pilon were very tunable, with error counts depending heavily on parameters, and their 'best' parameters were much better than their defaults.
* HyPo was not that tunable, with 'best' parameters barely better than defaults.
* The other three tools (NextPolish, Polypolish and Pypolca) were somewhat tunable - their 'best' parameters did better than defaults but not drastically so.
* Pypolca and Polypolish were still the best low-depth polishers after tuning, though Pilon and FMLRC2 did nearly as well.
* This parameter sweep is almost certainly overfitting to our data: _C. lari_, 18 errors, 5-20x depth.
* This sweep wasn't exhaustive. Especially for the higher-parameter tools (e.g. NextPolish) there could be even better parameter sets that I missed.
