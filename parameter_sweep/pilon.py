#!/usr/bin/env python3
"""
Pilon parameter sweep

parameters adjusted:
  --fix fixlist    A comma-separated list of categories of issues to try to fix: "snps": try to fix individual base errors; "indels": try to fix small indels; "gaps": try to fill gaps; "local": try to detect and fix local misassemblies; "all": all of the above (default); "bases": shorthand for "snps" and "indels" (for back compatibility); "none": none of the above; new fasta file will not be written. The following are experimental fix types: "amb": fix ambiguous bases in fasta output (to most likely alternative); "breaks": allow local reassembly to open new gaps (with "local"); "circles": try to close circlar elements when used with long corrected reads; "novel": assemble novel sequence from unaligned non-jump reads.
  --flank nbases    Controls how much of the well-aligned reads will be used; this many bases at each end of the good reads will be ignored (default 10).
  --K    Kmer size used by internal assembler (default 47).
  --mindepth depth    Variants (snps and indels) will only be called if there is coverage of good pairs at this depth or more; if this value is >= 1, it is an absolute depth, if it is a fraction < 1, then minimum depth is computed by multiplying this value by the mean coverage for the region, with a minumum value of 5 (default 0.1: min depth to call  is 10% of mean coverage or 5, whichever is greater).
  --minmq    Minimum alignment mapping quality for a read to count in pileups (default 0)
  --minqual    Minimum base quality to consider for pileups (default 0)
  --nostrays    Skip making a pass through the input BAM files to identify stray pairs, that is, those pairs in which both reads are aligned but not marked valid because they have inconsistent orientation or separation. Identifying stray pairs can help fill gaps and assemble larger insertions, especially of repeat content.  However, doing so sometimes consumes considerable memory.

parameters NOT adjusted:
  --variant Sets up heuristics for variant calling, as opposed to assembly improvement; equivalent to "--vcf --fix all,breaks".
  --chunksize Input FASTA elements larger than this will be processed in smaller pieces not to exceed this size (default 10000000).
  --diploid    Sample is from diploid organism; will eventually affect calling of heterozygous SNPs
  --dumpreads    Dump reads for local re-assemblies.
  --duplicates    Use reads marked as duplicates in the input BAMs (ignored by default).
  --iupac    Output IUPAC ambiguous base codes in the output FASTA file when appropriate.
  --nonpf    Use reads which failed sequencer quality filtering (ignored by default).
  --targets targetlist    Only process the specified target(s).  Targets are comma-separated, and each target is a fasta element name optionally followed by a base range. Example: "scaffold00001,scaffold00002:10000-20000" would result in processing all of scaffold00001 and coordinates 10000-20000 of scaffold00002. If "targetlist" is the name of a file, each line will be treated as a target specification.
  --verbose    More verbose output.
  --debug    Debugging output (implies verbose).
  --defaultqual qual    Assumes bases are of this quality if quals are no present in input BAMs (default 10).
  --gapmargin    Closed gaps must be within this number of bases of true size to be closed (100000)
  --mingap    Minimum size for unclosed gaps (default 10)

Skipping --defaultqual because my reads all have qscores. Skipping --gapmargin or --mapgap because
there are no gaps in these assemblies to be closed.
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
        errors = run_pilon_each_depth(parameters_str, description)
        if errors <= best_errors:
            best_errors = errors
            BEST_PARAMETERS = parameters
        count += 1
        if count == 1000:
            break


def run_pilon_each_depth(parameters_str, description):
    errors_05 = run_pilon(parameters_str, '05.0')
    errors_10 = run_pilon(parameters_str, '10.0')
    errors_15 = run_pilon(parameters_str, '15.0')
    errors_20 = run_pilon(parameters_str, '20.0')
    total_errors = errors_05 + errors_10 + errors_15 + errors_20
    with open(OUTPUT, 'at') as f:
        f.write(f'Pilon\t{parameters_str}\t{description}\t{errors_05}\t{errors_10}\t{errors_15}\t{errors_20}\t{total_errors}\n')
    return total_errors


def run_pilon(parameters, depth):
    delete_file('pilon.fasta')

    command = f'pilon {parameters} --genome {DRAFT} --frags alignments_{depth}.bam --output pilon'
    print(command)
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError:
        return float('inf')

    command = f'seqtk seq -U pilon.fasta > temp.fasta && mv temp.fasta pilon.fasta'
    print(command)
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError:
        return float('inf')

    if not file_exists('pilon.fasta'):
        return float('inf')

    command = f'compare_assemblies.py --aligner edlib {REF} pilon.fasta | grep -o "*" | wc -l'
    print(command)
    try:
        errors = subprocess.check_output(command, shell=True, text=True).strip()
    except subprocess.CalledProcessError:
        return float('inf')

    delete_file('pilon.fasta')
    return int(errors)


def align_reads(depth):
    reads_1 = BASE_DIR + depth + '/ATCC_35221_Campylobacter_lari/reads_1.fastq.gz'
    reads_2 = BASE_DIR + depth + '/ATCC_35221_Campylobacter_lari/reads_2.fastq.gz'

    command = f'bwa mem -t 24 {DRAFT} {reads_1} {reads_2} | samtools sort > alignments_{depth}.bam; samtools index alignments_{depth}.bam'
    print(command)
    subprocess.run(command, shell=True, check=True)


def get_parameters():
    yield ['--fix', 'all', '--flank', '10', '--K', '47', '--mindepth', '0.1', '--minmq', '0', '--minqual', '0'], 'defaults'
    yield ['--fix', 'bases', '--flank', '10', '--K', '47', '--mindepth', '0.1', '--minmq', '0', '--minqual', '0'], 'bases'
    while True:
        yield get_random_parameters(), 'random'
        yield mutate_parameters(), 'mutated'


def get_random_parameters():
    parameters = []
    if random.random() < 0.5:
         parameters += ['--fix', 'all']
    else:
         parameters += ['--fix', 'bases']
    parameters += ['--flank', str(random.randint(0, 20))]  # 0 to 20
    parameters += ['--K', str(random.randint(5, 49) * 2 + 1)]  # 11 to 99 (odd numbers only)
    if random.random() < 0.5:
         parameters += ['--mindepth', str(round(random.uniform(0, 1), 2))]  # 0 to 1
    else:
         parameters += ['--mindepth', str(round(random.uniform(1, 3) ** 3))]  # 1 to 27, weighted to the low end
    parameters += ['--minmq', str(random.randint(0, 60))]  # 0 to 60
    parameters += ['--minqual', str(random.randint(0, 40))]  # 0 to 40
    if random.random() < 0.5:
         parameters.append('--nostrays')
    return parameters


def mutate_parameters():
    global BEST_PARAMETERS
    parameters = list(BEST_PARAMETERS)

    # mutate --fix
    fix_index = parameters.index('--fix') + 1
    fix_val = parameters[fix_index]
    if random.random() < 0.25:
        fix_val = 'all' if random.random() < 0.5 else 'bases'

    # mutate --flank
    flank_index = parameters.index('--flank') + 1
    flank_val = int(parameters[flank_index])
    if random.random() < 0.25:
        flank_val = flank_val + random.randint(-3, 3)
        flank_val = max(flank_val, 0)

    # mutate --K
    k_index = parameters.index('--K') + 1
    k_val = int(parameters[k_index])
    if random.random() < 0.25:
        k_val = k_val + random.choice([-2, 2])
        k_val = max(k_val, 11)
        k_val = min(k_val, 99)

    # mutate --mindepth
    mindepth_index = parameters.index('--mindepth') + 1
    if '.' in parameters[mindepth_index]:
        mindepth_val = float(parameters[mindepth_index])
    else:
        mindepth_val = int(parameters[mindepth_index])
    if random.random() < 0.25:
        if isinstance(mindepth_val, float):
            mindepth_val = mindepth_val + random.uniform(-0.1, 0.1)
            mindepth_val = round(mindepth_val, 2)
            mindepth_val = max(mindepth_val, 0.0)
            mindepth_val = min(mindepth_val, 0.99)
        else:
            mindepth_val = int(mindepth_val) + random.choice([-1, 1])
            mindepth_val = max(mindepth_val, 1)

    # mutate --minmq
    minmq_index = parameters.index('--minmq') + 1
    minmq_val = int(parameters[minmq_index])
    if random.random() < 0.25:
        minmq_val = minmq_val + random.randint(-10, 10)
        minmq_val = max(minmq_val, 0)
        minmq_val = min(minmq_val, 60)

    # mutate --minqual
    minqual_index = parameters.index('--minqual') + 1
    minqual_val = int(parameters[minqual_index])
    if random.random() < 0.25:
        minqual_val = minqual_val + random.randint(-5, 5)
        minqual_val = max(minqual_val, 0)
        minqual_val = min(minqual_val, 40)

    return ['--fix', fix_val, '--flank', str(flank_val), '--K', str(k_val), '--mindepth', str(mindepth_val), '--minmq', str(minmq_val), '--minqual', str(minqual_val)]


def delete_file(filename):
    try:
        os.remove(filename)
    except OSError:
        pass


def file_exists(filename):
    return os.path.exists(filename) and os.path.getsize(filename) > 0


if __name__ == '__main__':
    main()
