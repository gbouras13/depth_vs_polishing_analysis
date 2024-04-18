#!/usr/bin/env python3
"""
This script produces analyses the errors from compare_assemblies.py run on a variety of depths (from 25.0 to 50.0) 
for the supplementary analysis added to the How Low Can You Go Manuscript.

I required all the assemblies and errors as described in the readme

It can be run like this to view the results directly in the terminal:
  analyse_errors.py -i assemblies -o errors_summary_output

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<https://www.gnu.org/licenses/>.
"""

import argparse
from argparse import RawTextHelpFormatter
import datetime
import gzip
import os
import re
import pathlib
import pandas as pd
import shutil
import subprocess
import sys
import textwrap


def get_input():
    usage = 'python3 analyse_errors.py ...'
    parser = argparse.ArgumentParser(description='parses viral nmpfs FASTA', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--indir', action="store", help='input directory with assemblies and errors',  required=True)
    parser.add_argument('-o', '--outdir', action="store", help='outdir',  required=True)
    parser.add_argument('-c', '--csv', action="store", help='ont_error_list.csv',  required=True)
    parser.add_argument('--min_depth', action="store",help='minimum depth (25.0x)', default=25.0)
    parser.add_argument('--max_depth', action="store",help='maximum depth (50.0x)',  default=50.0)
    parser.add_argument(
        "-f", "--force", help="Overwrites the output directory.", action="store_true"
    )

    args = parser.parse_args()

    return args


def main():

    args = get_input()


    if args.force == True:

        if os.path.isdir(args.outdir) == True:
            print(
                f"Removing output directory {args.outdir} as -f or --force was specified."
            )
            shutil.rmtree(args.outdir)
        elif os.path.isfile(args.outdir) == True:
            os.remove(args.outdir)
        else:
            print(
                f"--force was specified even though the output directory {args.outdir} does not already exist. Continuing."
            )
    else:
        if os.path.isdir(args.outdir) == True or os.path.isfile(args.outdir) == True:
            print(
                f"The output directory {args.outdir} already exists and force was not specified. Please specify -f or --force to overwrite it."
            )
            sys.exit()

    # make dir
    if os.path.isdir(args.outdir) == False:
        os.mkdir(args.outdir)


    ont_errors = pd.read_csv(args.csv, sep = ",")

    existing_errors = set(ont_errors['ref'].to_list())

    # Generate list of floats

    depth_list = [round(x, 1) for x in range(int(float(args.min_depth) * 10), int((float(args.max_depth)+0.1) * 10), 1)]
    # Convert each element to float
    depth_list = [float(x) / 10 for x in depth_list]

    genomes = ["ATCC_10708_Salmonella_enterica","ATCC_19119_Listeria_ivanovii","ATCC_35221_Campylobacter_lari", "ATCC_14035_Vibrio_cholerae", "ATCC_25922_Escherichia_coli", "ATCC_35897_Listeria_welshimeri", "ATCC_17802_Vibrio_parahaemolyticus", "ATCC_33560_Campylobacter_jejuni", "ATCC_BAA-679_Listeria_monocytogenes"]
    polishers = ["fmlrc2","hypo","nextpolish","pilon","polypolish-careful", "polypolish", "pypolca", "pypolca-careful"]

    data_list = []

    for depth in depth_list:

        for genome in genomes:

            for polisher in polishers:

                error_file = os.path.join(args.indir, str(depth), genome, f"{polisher}.errors")

                # Open the file and print each line
                with open(error_file, 'r') as file:
                    # for line in file:
                    #     print(line.strip())

                    lines = file.readlines()
                    
                    # the first line is the query -> so gap '-' in the query/first line = deletion
                    # the second line is the reference -> so gap '-' in the reference/second line = insertion
                    # third counts the errors
                    # fourth is empty

                    if len(lines) > 3:

                        for i in range(0, len(lines), 4):
                            batch = lines[i:i+4]
                            
                            if batch:
                                third_line = batch[2] 
                                count_errors = third_line.count('*')

                                first_line_parts = batch[0].split(':')
                                second_line_parts = batch[1].split(':')
                                if len(second_line_parts) > 1:
                                    chrom_region = second_line_parts[0].strip()
                                    insertion_errors = second_line_parts[1].count('-')
                                else:
                                    print("No colon found in line 2.")

                                if len(second_line_parts) > 1:
                                    deletion_errors = first_line_parts[1].count('-')
                                else:
                                    print("No colon found in line 1.")

                                error_type = "SNP"

                                snp_errors = count_errors - (deletion_errors + insertion_errors)

                                if deletion_errors == 0 and insertion_errors == 0:
                                    error_type = "SNP"
                                elif snp_errors == 0 and insertion_errors == 0:
                                    error_type = "Deletion"
                                elif snp_errors == 0 and deletion_errors == 0:
                                    error_type = "Insertion"
                                else: 
                                    error_type = "Mixed"

                                print('Error chrom region', chrom_region)
                                print('Error number', count_errors)
                                print('Error type', error_type)

                                ref = f"{genome} {chrom_region}"

                                if ref in existing_errors:
                                    existing_flag = True
                                else:
                                    existing_flag = False

                                # specifically if the error is a triple deletion/double -> smaller deletion, mark that as existing
                                # I have checked these manually too
                                if genome == "ATCC_35221_Campylobacter_lari" and chrom_region == "chromosome 548172-548202": # triple
                                    existing_flag = True
                                elif genome == "ATCC_33560_Campylobacter_jejuni" and chrom_region == "chromosome 148684-148714": #double
                                    existing_flag = True
                                elif  genome == "ATCC_19119_Listeria_ivanovii" and chrom_region == "chromosome 2671889-2671920": # 3 -> 2
                                    existing_flag = True
                                elif genome == "ATCC_19119_Listeria_ivanovii" and chrom_region == "chromosome 2671890-2671920": # 3 -> 1
                                    existing_flag = True


                                entry = {'genome': genome,
                                         'depth': depth,  
                                         'polisher': polisher,
                                         'chrom_region': chrom_region,
                                         'error_count': count_errors,
                                         'error_type':error_type,
                                         'existing': existing_flag}
                                         
                                data_list.append(entry)
                                

    df = pd.DataFrame(data_list)
    df.to_csv(os.path.join(args.outdir, "per_error_output.csv"), index=False)


    # count errors per polisher, depth

    # Group by 'polisher' and 'depth', then aggregate
    agg_df = df.groupby(['polisher', 'depth']).agg(
        total_SNP_errors=('error_count', lambda x: (x * (df['error_type'] == 'SNP')).sum()),
        total_Deletion_errors=('error_count', lambda x: (x * (df['error_type'] == 'Deletion')).sum()),
        total_Insertion_errors=('error_count', lambda x: (x * (df['error_type'] == 'Insertion')).sum()),
        total_Mixed_errors=('error_count', lambda x: (x * (df['error_type'] == 'Mixed')).sum()),
        Existing=('error_count', lambda x: (x * df['existing']).sum()),
        Introduced=('error_count', lambda x: (x * (~df['existing'])).sum())
    ).reset_index()

    agg_df.to_csv(os.path.join(args.outdir, "aggregated_output.csv"), index=False)


    



if __name__ == "__main__":
    main()


