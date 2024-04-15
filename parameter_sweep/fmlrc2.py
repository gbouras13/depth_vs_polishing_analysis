#!/usr/bin/env python3
"""
FMLRC2 parameter sweep

parameters adjusted:
    -k, --K <kmer_sizes>...                k-mer sizes for correction, can be specified multiple times (default: "-k 21 59")
    -m, --min_count <min_count>            absolute minimum k-mer count to consisder a path (default: 5)
    -f, --min_dynamic_count <min_frac>     dynamic minimum k-mer count fraction of median to consider a path (default: 0.1)

parameters NOT adjusted:
    -b, --begin_index <begin_id>           index of read to start with (default: 0)
    -B, --branch_factor <branch_factor>    branching factor for correction, scaled by k (default: 4.0)
    -C, --cache_size <cache_size>          the length of k-mer to precompute in cache (default: 8)
    -e, --end_index <end_id>               index of read to end with (default: end of file)
    -t, --threads <threads>                number of correction threads (default: 1)

I skipped -b and -e because they just limit which reads are used. I skipped -C because it seems to
be a CPU vs RAM performance thing. I skipped -t because it's a performance thing. And I while I
don't understand what -B does, it didn't change the results when I adjusted it, so I skipped it.
"""

import random
import subprocess


BASE_DIR = '/home/wickr/2024-01_low-depth_Illumina_polishing/'
REF = BASE_DIR + 'references/ATCC_35221_Campylobacter_lari.fasta'
DRAFT = BASE_DIR + 'drafts/ATCC_35221_Campylobacter_lari.fasta'
OUTPUT = 'results.tsv'

BEST_PARAMETERS = None


def main():
    global BEST_PARAMETERS
    for depth in ['05.0', '10.0', '15.0', '20.0']:
        build_bwt(depth)
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
        errors = run_fmlrc2_each_depth(parameters_str, description)
        if errors <= best_errors:
            best_errors = errors
            BEST_PARAMETERS = parameters
        count += 1
        if count == 1000:
            break


def run_fmlrc2_each_depth(parameters_str, description):
    errors_05 = run_fmlrc2(parameters_str, '05.0')
    errors_10 = run_fmlrc2(parameters_str, '10.0')
    errors_15 = run_fmlrc2(parameters_str, '15.0')
    errors_20 = run_fmlrc2(parameters_str, '20.0')
    total_errors = errors_05 + errors_10 + errors_15 + errors_20
    with open(OUTPUT, 'at') as f:
        f.write(f'FMLRC2\t{parameters_str}\t{description}\t{errors_05}\t{errors_10}\t{errors_15}\t{errors_20}\t{total_errors}\n')
    return total_errors


def run_fmlrc2(parameters, depth):
    command = f'fmlrc2 {parameters} -t 24 comp_msbwt_{depth}.npy {DRAFT} fmlrc2.fasta'
    print(command)
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError:
        return float('inf')
    command = f'compare_assemblies.py --aligner edlib {REF} fmlrc2.fasta | grep -o "*" | wc -l'
    print(command)
    try:
        errors = subprocess.check_output(command, shell=True, text=True).strip()
    except subprocess.CalledProcessError:
        return float('inf')
    return int(errors)


def build_bwt(depth):
    reads_1=BASE_DIR + depth + '/ATCC_35221_Campylobacter_lari/reads_1.fastq.gz'
    reads_2=BASE_DIR + depth + '/ATCC_35221_Campylobacter_lari/reads_2.fastq.gz'
    command = f"gunzip -c {reads_1} {reads_2} | awk 'NR % 4 == 2' | tr NT TN | ropebwt2 -LR | tr NT TN | fmlrc2-convert comp_msbwt_{depth}.npy"
    print(command)
    subprocess.run(command, shell=True, check=True)


def get_parameters():
    yield ['-k', '21', '59', '-m', '5', '-f', '0.1'], 'defaults'
    yield ['-k', '21', '59', '80', '-m', '0', '-f', '0.1'], 'eukaryote mode'
    while True:
        yield get_random_parameters(), 'random'
        yield mutate_parameters(), 'mutated'


def get_random_parameters():
    kmers = sorted([round(int(random.uniform(3, 10) ** 2))
                    for _ in range(random.randint(1, 3))])
    parameters = ['-k'] + [str(k) for k in kmers]
    parameters += ['-m', str(round(random.uniform(0, 4) ** 3))]  # 0 to 64, weighted to the low end
    parameters += ['-f', str(round(random.uniform(0, 1) ** 3, 2))]  # 0 to 1, weighted to the low end
    return parameters


def mutate_parameters():
    global BEST_PARAMETERS
    parameters = list(BEST_PARAMETERS)

    # mutate -m
    m_index = parameters.index('-m') + 1
    m_val = int(parameters[m_index])
    if random.random() < 0.25:
        m_val = m_val + random.choice([-1, 1])
        m_val = max(m_val, 0)

    # mutate -f
    f_index = parameters.index('-f') + 1
    f_val = float(parameters[f_index])
    if random.random() < 0.25:
        f_val = f_val + random.uniform(-0.1, 0.1)
        f_val = max(f_val, 0.0)

    # mutate -k
    k_vals = [int(k) for k in parameters[1:m_index-1]]
    if random.random() < 0.1:
        k_vals.append(random.randint(10, 100))
    if random.random() < 0.1 and len(k_vals) > 1:
        random.shuffle(k_vals)
        k_vals.pop()
    for i in range(len(k_vals)):
        if random.random() < 0.25:
            k_vals[i] += random.choice([-2, -1, 1, 2])
    k_vals = sorted(list(set(k_vals)))
    k_vals = [str(k) for k in k_vals]

    return ['-k'] + k_vals + ['-m', str(max(0, m_val)), '-f', str(round(max(0.0, f_val), 2))]


if __name__ == '__main__':
    main()
