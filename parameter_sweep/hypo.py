#!/usr/bin/env python3
"""
HyPo parameter sweep

parameters adjusted:
  -m, --match-sr <int>         Score for matching bases for short reads. [Default] 5.
  -x, --mismatch-sr <int>      Score for mismatching bases for short reads. [Default] -4.
  -g, --gap-sr <int>           Gap penalty for short reads (must be negative). [Default] -8.
  -n, --ned-th <int>           Threshold for Normalised Edit Distance of long arms allowed in a window (in %). Higher number means more arms allowed which may slow down the execution. [Default] 20.
  -q, --qual-map-th <int>      Threshold for mapping quality of reads. The reads with mapping quality below this threshold will not be taken into consideration. [Default] 2.

parameters NOT adjusted:
  -t, --threads <int>          Number of threads. [Default] 1.
  -p, --processing-size <int>  Number of contigs to be processed in one batch. Lower value means less memory usage but slower speed. [Default] All the contigs in the draft.
  -k, --kind-sr <str>          Kind of the short reads. [Valid Values] sr (Corresponding to NGS reads like Illumina reads) ccs (Corresponding to HiFi reads like PacBio CCS reads) [Default] sr.
  -M, --match-lr <int>         Score for matching bases for long reads. [Default] 3.
  -X, --mismatch-lr <int>      Score for mismatching bases for long reads. [Default] -5.
  -G, --gap-lr <int>           Gap penalty for long reads (must be negative). [Default] -4.
  -i, --intermed               Store or use (if already exist) the intermediate files. [Currently, only Solid kmers are stored as an intermediate file.].

I skipped -t, -p and -i because they just seem to be performance related. I skipped -k because I'm
only using Illumina reads. And I skipped -M, -X and -G because they are for long reads (which I'm
not using).
"""

import os
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
    for depth in ['05.0', '10.0', '15.0', '20.0']:
        align_reads(depth)
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
        errors = run_hypo_each_depth(parameters_str, description)
        if errors <= best_errors:
            best_errors = errors
            BEST_PARAMETERS = parameters
        count += 1
        if count == 1000:
            break


def run_hypo_each_depth(parameters_str, description):
    errors_05 = run_hypo(parameters_str, '05.0')
    errors_10 = run_hypo(parameters_str, '10.0')
    errors_15 = run_hypo(parameters_str, '15.0')
    errors_20 = run_hypo(parameters_str, '20.0')
    total_errors = errors_05 + errors_10 + errors_15 + errors_20
    with open(OUTPUT, 'at') as f:
        f.write(f'HyPo\t{parameters_str}\t{description}\t{errors_05}\t{errors_10}\t{errors_15}\t{errors_20}\t{total_errors}\n')
    return total_errors


def run_hypo(parameters, depth):
    try:
        os.remove('hypo.fasta')
    except OSError:
        pass
    integer_depth = int(float(depth))
    command = f'hypo {parameters} -d {DRAFT} -r @read_filenames_{depth}.txt -s 1513368 -c {integer_depth} -b alignments_{depth}.bam -t 24 -o hypo.fasta'
    print(command)
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError:
        shutil.rmtree('aux', ignore_errors=True)
    shutil.rmtree('aux', ignore_errors=True)
    if not os.path.exists('hypo.fasta'):
        return float('inf')
    command = f'compare_assemblies.py --aligner edlib {REF} hypo.fasta | grep -o "*" | wc -l'
    print(command)
    try:
        errors = subprocess.check_output(command, shell=True, text=True).strip()
    except subprocess.CalledProcessError:
        return float('inf')
    return int(errors)


def align_reads(depth):
    reads_1 = BASE_DIR + depth + '/ATCC_35221_Campylobacter_lari/reads_1.fastq.gz'
    reads_2 = BASE_DIR + depth + '/ATCC_35221_Campylobacter_lari/reads_2.fastq.gz'

    command = f'bwa mem -t 24 {DRAFT} {reads_1} {reads_2} | samtools sort > alignments_{depth}.bam; samtools index alignments_{depth}.bam'
    print(command)
    subprocess.run(command, shell=True, check=True)

    command = f'echo -e "{reads_1}\n{reads_2}" > read_filenames_{depth}.txt'
    print(command)
    subprocess.run(command, shell=True, check=True)


def get_parameters():
    yield ['-m', '5', '-x', '-4', '-g', '-8', '-n', '20', '-q', '2'], 'defaults'
    while True:
        yield get_random_parameters(), 'random'
        yield mutate_parameters(), 'mutated'


def get_random_parameters():
    parameters = ['-m', str(random.randint(0, 10))]  # 0 to 10
    parameters += ['-x', str(random.randint(-10, 0))]  # -10 to 0
    parameters += ['-g', str(random.randint(-10, 0))]  # -10 to 0
    parameters += ['-n', str(round(random.uniform(0, 10) ** 2))]  # 0 to 100, weighted to the low end
    parameters += ['-q', str(random.randint(0, 60))]  # 0 to 60
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

    # mutate -x
    x_index = parameters.index('-x') + 1
    x_val = int(parameters[x_index])
    if random.random() < 0.25:
        x_val = x_val + random.choice([-1, 1])
        x_val = min(x_val, 0)

    # mutate -g
    g_index = parameters.index('-g') + 1
    g_val = int(parameters[g_index])
    if random.random() < 0.25:
        g_val = g_val + random.choice([-1, 1])
        g_val = min(g_val, 0)

    # mutate -n
    n_index = parameters.index('-n') + 1
    n_val = int(parameters[n_index])
    if random.random() < 0.25:
        n_val = n_val + random.randint(-10, 10)
        n_val = max(n_val, 0)
        n_val = min(n_val, 100)

    # mutate -q
    q_index = parameters.index('-q') + 1
    q_val = int(parameters[q_index])
    if random.random() < 0.25:
        q_val = q_val + random.randint(-10, 10)
        q_val = max(q_val, 0)
        q_val = min(q_val, 60)

    return ['-m', str(m_val), '-x', str(x_val), '-g', str(g_val), '-n', str(n_val), '-q', str(q_val)]


if __name__ == '__main__':
    main()
