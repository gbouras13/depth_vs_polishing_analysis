#!/usr/bin/env python3
"""
Polypolish parameter sweep

parameters adjusted:
  -i, --fraction_invalid <FRACTION_INVALID>  A base must make up less than this fraction of the read depth to be considered invalid [default: 0.2]
  -v, --fraction_valid <FRACTION_VALID>      A base must make up at least this fraction of the read depth to be considered valid [default: 0.5]
  -m, --max_errors <MAX_ERRORS>              Ignore alignments with more than this many mismatches and indels [default: 10]
  -d, --min_depth <MIN_DEPTH>                A base must occur at least this many times in the pileup to be considered valid [default: 5]
      --careful                              Ignore any reads with multiple alignments

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
        align_and_filter(depth)
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
        errors = run_polypolish_each_depth(parameters_str, description)
        if errors <= best_errors:
            best_errors = errors
            BEST_PARAMETERS = parameters
        count += 1
        if count == 1000:
            break


def run_polypolish_each_depth(parameters_str, description):
    errors_05 = run_polypolish(parameters_str, '05.0')
    errors_10 = run_polypolish(parameters_str, '10.0')
    errors_15 = run_polypolish(parameters_str, '15.0')
    errors_20 = run_polypolish(parameters_str, '20.0')
    total_errors = errors_05 + errors_10 + errors_15 + errors_20
    with open(OUTPUT, 'at') as f:
        f.write(f'Polypolish\t{parameters_str}\t{description}\t{errors_05}\t{errors_10}\t{errors_15}\t{errors_20}\t{total_errors}\n')
    return total_errors


def run_polypolish(parameters, depth):
    command = f'polypolish polish {parameters} {DRAFT} filtered_1_{depth}.sam filtered_2_{depth}.sam 1> polypolish.fasta'
    print(command)
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError:
        return float('inf')
    command = f'compare_assemblies.py --aligner edlib {REF} polypolish.fasta | grep -o "*" | wc -l'
    print(command)
    try:
        errors = subprocess.check_output(command, shell=True, text=True).strip()
    except subprocess.CalledProcessError:
        return float('inf')
    return int(errors)


def align_and_filter(depth):
    reads_1 = BASE_DIR + depth + '/ATCC_35221_Campylobacter_lari/reads_1.fastq.gz'
    reads_2 = BASE_DIR + depth + '/ATCC_35221_Campylobacter_lari/reads_2.fastq.gz'

    command = f'bwa mem -t 24 -a {DRAFT} {reads_1} > alignments_1.sam'
    print(command)
    subprocess.run(command, shell=True, check=True)
    
    command = f'bwa mem -t 24 -a {DRAFT} {reads_2} > alignments_2.sam'
    print(command)
    subprocess.run(command, shell=True, check=True)
    
    command = f'polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1_{depth}.sam --out2 filtered_2_{depth}.sam'
    print(command)
    subprocess.run(command, shell=True, check=True)


def get_parameters():
    yield ['-i', '0.2', '-v', '0.5', '-m', '10', '-d', '5'], 'defaults'
    yield ['-i', '0.2', '-v', '0.5', '-m', '10', '-d', '5', '--careful'], 'careful mode'
    while True:
        yield get_random_parameters(), 'random'
        yield mutate_parameters(), 'mutated'


def get_random_parameters():
    i_and_v = sorted([round(random.uniform(0, 0.8), 2), round(random.uniform(0.2, 1.0), 2)])
    parameters = ['-i', str(i_and_v[0])]
    parameters += ['-v', str(i_and_v[1])]
    parameters += ['-m', str(random.randint(0, 20))]
    parameters += ['-d', str(random.randint(1, 20))]
    if random.random() < 0.5:
         parameters.append('--careful')
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

    # mutate -d
    d_index = parameters.index('-d') + 1
    d_val = int(parameters[d_index])
    if random.random() < 0.25:
        d_val = d_val + random.choice([-1, 1])
        d_val = max(d_val, 1)

    # mutate -i and -v
    i_index = parameters.index('-i') + 1
    i_val = float(parameters[i_index])
    if random.random() < 0.25:
        i_val = i_val + random.uniform(-0.1, 0.1)
        i_val = max(i_val, 0.0)
        i_val = min(i_val, 0.8)
    v_index = parameters.index('-v') + 1
    v_val = float(parameters[v_index])
    if random.random() < 0.25:
        v_val = v_val + random.uniform(-0.1, 0.1)
        v_val = max(v_val, 0.2)
        v_val = min(v_val, 1.0)
    i_val, v_val = sorted([i_val, v_val])

    parameters = ['-i', str(round(i_val, 2)), '-v', str(round(v_val, 2)), '-m', str(round(m_val, 2)), '-d', str(round(d_val, 2))]
    if random.random() < 0.5:
         parameters.append('--careful')

    return parameters


if __name__ == '__main__':
    main()
