#!/usr/bin/env python3
"""
Pypolca parameter sweep

parameters adjusted:
  --min_alt INTEGER        Minimum alt allele count to make a change [default: 2]
  --min_ratio FLOAT        Minimum alt allele to ref allele ratio to make a change  [default: 2.0]
  --careful                Equivalent to --min_alt 4 --min_ratio 3

parameters NOT adjusted:
  -t, --threads INTEGER    Number of threads.  [default: 1]
  -o, --output PATH        Output directory path  [default: output_polca]
  -f, --force              Force overwrites the output directory
  -n, --no_polish          do not polish, just create vcf file, evaluate the assembly and exit
  -m, --memory_limit TEXT  Memory per thread to use in samtools sort, set to 2G or more for large genomes  [default: 2G]
  -p, --prefix TEXT        prefix for output files  [default: polca]

I skipped -b and -e because they just limit which reads are used. I skipped -C because it seems to
be a CPU vs RAM performance thing. I skipped -t because it's a performance thing. And I while I
don't understand what -B does, it didn't change the results when I adjusted it, so I skipped it.
"""

import random
import shutil
import subprocess


BASE_DIR = '/home/wickr/2024-01_low-depth_Illumina_polishing/'
REF = BASE_DIR + 'references/ATCC_35221_Campylobacter_lari.fasta'
DRAFT = BASE_DIR + 'drafts/ATCC_35221_Campylobacter_lari.fasta'
OUTPUT = 'results.tsv'

BEST_PARAMETERS = None


def main():
    global BEST_PARAMETERS
    with open(OUTPUT, 'at') as f:
        f.write(f'tool\tparameters\tdescription\terrors_05\terrors_10\terrors_15\terrors_20\ttotal_errors\n')
    tried_parameters = set()
    best_errors = float('inf')
    count = 0
    for parameters, description in get_parameters():
        parameters_str = ' '.join(parameters)
        if parameters_str in tried_parameters:
            continue
        tried_parameters.add(parameters_str)
        errors = run_pypolca_each_depth(parameters_str, description)
        if errors <= best_errors:
            best_errors = errors
            BEST_PARAMETERS = parameters
        count += 1
        if count == 1000:
            break


def run_pypolca_each_depth(parameters_str, description):
    errors_05 = run_pypolca(parameters_str, '05.0')
    errors_10 = run_pypolca(parameters_str, '10.0')
    errors_15 = run_pypolca(parameters_str, '15.0')
    errors_20 = run_pypolca(parameters_str, '20.0')
    total_errors = errors_05 + errors_10 + errors_15 + errors_20
    with open(OUTPUT, 'at') as f:
        f.write(f'Pypolca\t{parameters_str}\t{description}\t{errors_05}\t{errors_10}\t{errors_15}\t{errors_20}\t{total_errors}\n')
    return total_errors


def run_pypolca(parameters, depth):
    reads_1=BASE_DIR + depth + '/ATCC_35221_Campylobacter_lari/reads_1.fastq.gz'
    reads_2=BASE_DIR + depth + '/ATCC_35221_Campylobacter_lari/reads_2.fastq.gz'
    command = f'pypolca run {parameters} -a {DRAFT} -1 {reads_1} -2 {reads_2} -t 24 -o pypolca_temp'
    print(command)
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError:
        return float('inf')
    command = f'compare_assemblies.py --aligner edlib {REF} pypolca_temp/polca_corrected.fasta | grep -o "*" | wc -l'
    print(command)
    try:
        errors = subprocess.check_output(command, shell=True, text=True).strip()
    except subprocess.CalledProcessError:
        return float('inf')
    shutil.rmtree('pypolca_temp')
    return int(errors)


def get_parameters():
    yield ['--min_alt', '2', '--min_ratio', '2.0'], 'defaults'
    yield ['--min_alt', '4', '--min_ratio', '3.0'], 'careful mode'
    while True:
        yield get_random_parameters(), 'random'
        yield mutate_parameters(), 'mutated'


def get_random_parameters():
    parameters = ['--min_alt', str(round(random.uniform(1, 3) ** 3))]  # 1 to 27, weighted to the low end
    parameters += ['--min_ratio', str(round(random.uniform(0.8, 3.0) ** 3, 2))]  # 0.51 to 27.0, weighted to the low end
    return parameters


def mutate_parameters():
    global BEST_PARAMETERS
    parameters = list(BEST_PARAMETERS)

    # mutate --min_alt
    min_alt_index = parameters.index('--min_alt') + 1
    min_alt_val = int(parameters[min_alt_index])
    if random.random() < 0.5:
        min_alt_val = min_alt_val + random.choice([-1, 1])
        min_alt_val = max(min_alt_val, 1)

    # mutate --min_ratio
    min_ratio_index = parameters.index('--min_ratio') + 1
    min_ratio_val = float(parameters[min_ratio_index])
    if random.random() < 0.5:
        min_ratio_val = min_ratio_val + random.uniform(-0.1, 0.1)
        min_ratio_val = max(min_ratio_val, 0.0)

    return ['--min_alt', str(min_alt_val), '--min_ratio', str(round(min_ratio_val, 2))]


if __name__ == '__main__':
    main()
