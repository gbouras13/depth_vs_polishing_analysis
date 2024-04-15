#!/usr/bin/env python3
"""
NextPolish parameter sweep

parameters adjusted:
  -min_map_quality N           skip the mapped read with mapping quality < N. (default: 0)
  -max_ins_fold_sgs N          skip the paired-end read with insert size > N * estimated average insert size. (default: 5)
  -max_clip_ratio_sgs F        skip the mapped read with clipped length > F * full length, used for bam_sgs. (default: 0.15)
  -trim_len_edge N             trimed length at the two edges of a alignment. (default: 2)
  -ext_len_edge N              extened length at the two edges of a low quality region. (default: 2)
  -indel_balance_factor_sgs F  a factor to control the ratio between indels, larger factor will produced more deletions, and vice versa. (default: 0.5)
  -min_count_ratio_skip F      skip a site if the fraction of the most genotype > F. (default: 0.8)
  -max_len_kmer N              maximum length requirement of a polished kmer, longer kmers will be splited. (default: 50)
  -min_len_inter_kmer N        minimum interval length between two adjacent kmers, shorter interval length will be merged. (default: 5)
  -max_count_kmer N            read up to this count of observed kmers for a polished kmer. (default: 50)

parameters NOT adjusted:
  -p N, --process N            number of processes used for polishing. (default: 10)
  -count_read_ins_sgs N        read N reads to estimate the insert size of paired-end reads. (default: 10000)
  -max_ins_len_sgs N           skip the paired-end read with insert size > N. (default: 10000)
  -max_clip_ratio_lgs F        skip the mapped read with clipped length > F * full length, used for bam_lgs. (default: 0.4)
  -min_len_ldr N               minimum length requirement of a low depth region, which will be further processed using bam_lgs. (default: 3)
  -ploidy N                    set the ploidy of the sample of this genome. (default: 2)
  -max_variant_count_lgs N     exclude long reads with more than N variable sites, it is approximately equivalent to total error bases in the long read. (default: 150k)
  -indel_balance_factor_lgs F  a factor to control the ratio between indels, larger factor will produced more deletions, and vice versa. (default: 0.33)
  -min_depth_snp N             recall snps using bam_lgs if the total depth of this site in bam_sgs < N. (default: 3)
  -min_count_snp N             recall snps using bam_lgs if the count of this snp in bam_sgs < N. (default: 5)
  -min_count_snp_link N        find a snp linkage using bam_lgs if the count of this linkage in bam_sgs < N. (default: 5)
  -max_indel_factor_lgs F       recall indels with bam_sgs if the count of the second most genotype > F * the count of the most genotype when the most genotype is different with ref in bam_lgs. (default: 0.21)
  -max_snp_factor_lgs F        recall snps with bam_lgs if the count of the second most genotype > F * the count of the most genotype when the most genotype is different with ref. (default: 0.53)
  -min_snp_factor_sgs F        skip a snp if the count of the second most genotype < F * the count of the most genotype. (default: 0.34)

Skipped -max_ins_len_sgs because it seems redudant with -max_ins_fold_sgs. Skipped
-max_clip_ratio_lgs and min_len_ldr because they're for long reads. Skipped -ploidy,
-max_variant_count_lgs, -indel_balance_factor_lgs, -min_depth_snp, -min_count_snp,
-min_count_snp_link, -max_indel_factor_lgs, -max_snp_factor_lgs and -min_snp_factor_sgs
because they are for a task I'm not running.

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
        errors = run_nextpolish_each_depth(parameters_str, description)
        if errors <= best_errors:
            best_errors = errors
            BEST_PARAMETERS = parameters
        count += 1
        if count == 1000:
            break


def run_nextpolish_each_depth(parameters_str, description):
    errors_05 = run_nextpolish(parameters_str, '05.0')
    errors_10 = run_nextpolish(parameters_str, '10.0')
    errors_15 = run_nextpolish(parameters_str, '15.0')
    errors_20 = run_nextpolish(parameters_str, '20.0')
    total_errors = errors_05 + errors_10 + errors_15 + errors_20
    with open(OUTPUT, 'at') as f:
        f.write(f'NextPolish\t{parameters_str}\t{description}\t{errors_05}\t{errors_10}\t{errors_15}\t{errors_20}\t{total_errors}\n')
    return total_errors


def run_nextpolish(parameters, depth):
    reads_1 = BASE_DIR + depth + '/ATCC_35221_Campylobacter_lari/reads_1.fastq.gz'
    reads_2 = BASE_DIR + depth + '/ATCC_35221_Campylobacter_lari/reads_2.fastq.gz'

    delete_file('nextpolish.fasta')

    command = f'python /home/wickr/programs/NextPolish/lib/nextpolish1.py {parameters} -g {DRAFT} -t 1 -p 24 -s sgs.sort_{depth}.bam > nextpolish_temp.fasta'
    print(command)
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError:
        return float('inf')

    command = f'bwa index nextpolish_temp.fasta && bwa mem -t 24 nextpolish_temp.fasta {reads_1} {reads_2} | samtools view --threads 3 -F 0x4 -b - | samtools fixmate -m --threads 3  - - | samtools sort -m 2g --threads 5 - | samtools markdup --threads 5 -r - sgs.sort_temp.bam && samtools index sgs.sort_temp.bam && samtools faidx nextpolish_temp.fasta'
    print(command)
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError:
        return float('inf')

    command = f'python /home/wickr/programs/NextPolish/lib/nextpolish1.py -g nextpolish_temp.fasta -t 2 -p 24 -s sgs.sort_temp.bam | seqkit sort -l -r | seqtk seq -U > nextpolish.fasta'
    print(command)
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError:
        return float('inf')

    if not file_exists('nextpolish.fasta'):
        return float('inf')

    delete_file('nextpolish_temp.fasta')
    delete_file('nextpolish_temp.fasta.fai')
    delete_file('nextpolish_temp.fasta.amb')
    delete_file('nextpolish_temp.fasta.ann')
    delete_file('nextpolish_temp.fasta.bwt')
    delete_file('nextpolish_temp.fasta.pac')
    delete_file('nextpolish_temp.fasta.sa')
    delete_file('sgs.sort_temp.bam')
    delete_file('sgs.sort_temp.bam.bai')
    
    command = f'compare_assemblies.py --aligner edlib {REF} nextpolish.fasta | grep -o "*" | wc -l'
    print(command)
    try:
        errors = subprocess.check_output(command, shell=True, text=True).strip()
    except subprocess.CalledProcessError:
        return float('inf')

    delete_file('nextpolish.fasta')
    return int(errors)


def align_reads(depth):
    reads_1 = BASE_DIR + depth + '/ATCC_35221_Campylobacter_lari/reads_1.fastq.gz'
    reads_2 = BASE_DIR + depth + '/ATCC_35221_Campylobacter_lari/reads_2.fastq.gz'

    command = f'bwa mem -t 24 {DRAFT} {reads_1} {reads_2} | samtools view --threads 3 -F 0x4 -b - | samtools fixmate -m --threads 3  - - | samtools sort -m 2g --threads 5 - | samtools markdup --threads 5 -r - sgs.sort_{depth}.bam && samtools index sgs.sort_{depth}.bam'
    print(command)
    subprocess.run(command, shell=True, check=True)


def get_parameters():
    yield ['-min_map_quality', '0', '-max_ins_fold_sgs', '5', '-max_clip_ratio_sgs', '0.15', '-trim_len_edge', '2',
           '-ext_len_edge', '2', '-indel_balance_factor_sgs', '0.5', '-min_count_ratio_skip', '0.8',
           '-max_len_kmer', '50', '-min_len_inter_kmer', '5', '-max_count_kmer', '50'], 'defaults'
    while True:
        yield get_random_parameters(), 'random'
        yield mutate_parameters(), 'mutated'


def get_random_parameters():
    parameters = ['-min_map_quality', str(random.randint(0, 60))]  # 0 to 60
    parameters += ['-max_ins_fold_sgs', str(round(random.uniform(1, 4) ** 2))]  # 1 to 16, weighted to the low end
    parameters += ['-max_clip_ratio_sgs', str(round(random.uniform(0, 1) ** 3, 2))]  # 0 to 1, weighted to the low end
    parameters += ['-trim_len_edge', str(round(random.uniform(0, 3) ** 3))]  # 0 to 27, weighted to the low end
    parameters += ['-ext_len_edge', str(round(random.uniform(0, 3) ** 3))]  # 0 to 27, weighted to the low end
    parameters += ['-indel_balance_factor_sgs', str(round(random.uniform(0, 1), 2))]  # 0 to 1
    parameters += ['-min_count_ratio_skip', str(round(random.uniform(0, 1) ** 0.5, 2))]  # 0 to 1, weighted to the high end
    parameters += ['-max_len_kmer', str(random.randint(10, 100))]  # 10 to 100
    parameters += ['-min_len_inter_kmer', str(random.randint(1, 10))]  # 1 to 10
    parameters += ['-max_count_kmer', str(random.randint(1, 100))]  # 1 to 100
    return parameters


def mutate_parameters():
    global BEST_PARAMETERS
    parameters = list(BEST_PARAMETERS)

    # mutate -min_map_quality
    min_map_quality_index = parameters.index('-min_map_quality') + 1
    min_map_quality_val = int(parameters[min_map_quality_index])
    if random.random() < 0.15:
        min_map_quality_val = min_map_quality_val + random.randint(-10, 10)
        min_map_quality_val = max(min_map_quality_val, 0)
        min_map_quality_val = min(min_map_quality_val, 60)

    # mutate -max_ins_fold_sgs
    max_ins_fold_sgs_index = parameters.index('-max_ins_fold_sgs') + 1
    max_ins_fold_sgs_val = int(parameters[max_ins_fold_sgs_index])
    if random.random() < 0.15:
        max_ins_fold_sgs_val = max_ins_fold_sgs_val + random.randint(-2, 2)
        max_ins_fold_sgs_val = max(max_ins_fold_sgs_val, 1)

    # mutate -max_clip_ratio_sgs
    max_clip_ratio_sgs_index = parameters.index('-max_clip_ratio_sgs') + 1
    max_clip_ratio_sgs_val = float(parameters[max_clip_ratio_sgs_index])
    if random.random() < 0.15:
        max_clip_ratio_sgs_val = max_clip_ratio_sgs_val + random.uniform(-0.1, 0.1)
        max_clip_ratio_sgs_val = round(max_clip_ratio_sgs_val, 2)
        max_clip_ratio_sgs_val = max(max_clip_ratio_sgs_val, 0.0)
        max_clip_ratio_sgs_val = min(max_clip_ratio_sgs_val, 1.0)

    # mutate -trim_len_edge
    trim_len_edge_index = parameters.index('-trim_len_edge') + 1
    trim_len_edge_val = int(parameters[trim_len_edge_index])
    if random.random() < 0.15:
        trim_len_edge_val = trim_len_edge_val + random.randint(-2, 2)
        trim_len_edge_val = max(trim_len_edge_val, 0)

    # mutate -ext_len_edge
    ext_len_edge_index = parameters.index('-ext_len_edge') + 1
    ext_len_edge_val = int(parameters[ext_len_edge_index])
    if random.random() < 0.15:
        ext_len_edge_val = ext_len_edge_val + random.randint(-2, 2)
        ext_len_edge_val = max(ext_len_edge_val, 0)

    # mutate -indel_balance_factor_sgs
    indel_balance_factor_sgs_index = parameters.index('-indel_balance_factor_sgs') + 1
    indel_balance_factor_sgs_val = float(parameters[indel_balance_factor_sgs_index])
    if random.random() < 0.15:
        indel_balance_factor_sgs_val = indel_balance_factor_sgs_val + random.uniform(-0.1, 0.1)
        indel_balance_factor_sgs_val = round(indel_balance_factor_sgs_val, 2)
        indel_balance_factor_sgs_val = max(indel_balance_factor_sgs_val, 0.0)
        indel_balance_factor_sgs_val = min(indel_balance_factor_sgs_val, 1.0)

    # mutate -min_count_ratio_skip
    min_count_ratio_skip_index = parameters.index('-min_count_ratio_skip') + 1
    min_count_ratio_skip_val = float(parameters[min_count_ratio_skip_index])
    if random.random() < 0.15:
        min_count_ratio_skip_val = min_count_ratio_skip_val + random.uniform(-0.1, 0.1)
        min_count_ratio_skip_val = round(min_count_ratio_skip_val, 2)
        min_count_ratio_skip_val = max(min_count_ratio_skip_val, 0.0)
        min_count_ratio_skip_val = min(min_count_ratio_skip_val, 1.0)

    # mutate -max_len_kmer
    max_len_kmer_index = parameters.index('-max_len_kmer') + 1
    max_len_kmer_val = int(parameters[max_len_kmer_index])
    if random.random() < 0.15:
        max_len_kmer_val = max_len_kmer_val + random.randint(-10, 10)
        max_len_kmer_val = max(max_len_kmer_val, 10)

    # mutate -min_len_inter_kmer
    min_len_inter_kmer_index = parameters.index('-min_len_inter_kmer') + 1
    min_len_inter_kmer_val = int(parameters[min_len_inter_kmer_index])
    if random.random() < 0.15:
        min_len_inter_kmer_val = min_len_inter_kmer_val + random.randint(-2, 2)
        min_len_inter_kmer_val = max(min_len_inter_kmer_val, 1)

    # mutate -max_count_kmer
    max_count_kmer_index = parameters.index('-max_count_kmer') + 1
    max_count_kmer_val = int(parameters[max_count_kmer_index])
    if random.random() < 0.15:
        max_count_kmer_val = max_count_kmer_val + random.randint(-10, 10)
        max_count_kmer_val = max(max_count_kmer_val, 1)

    return ['-min_map_quality', str(min_map_quality_val),
            '-max_ins_fold_sgs', str(max_ins_fold_sgs_val),
            '-max_clip_ratio_sgs', str(max_clip_ratio_sgs_val),
            '-trim_len_edge', str(trim_len_edge_val),
            '-ext_len_edge', str(ext_len_edge_val),
            '-indel_balance_factor_sgs', str(indel_balance_factor_sgs_val),
            '-min_count_ratio_skip', str(min_count_ratio_skip_val),
            '-max_len_kmer', str(max_len_kmer_val),
            '-min_len_inter_kmer', str(min_len_inter_kmer_val),
            '-max_count_kmer', str(max_count_kmer_val)]


def delete_file(filename):
    try:
        os.remove(filename)
    except OSError:
        pass


def file_exists(filename):
    return os.path.exists(filename) and os.path.getsize(filename) > 0


if __name__ == '__main__':
    main()
