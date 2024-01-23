# Use the compare_assemblies.py script to produce FASTA files of the error
# locations in the reference sequences:
rm -f *_error_locations.fasta
for a in ATCC_10708_Salmonella_enterica ATCC_14035_Vibrio_cholerae ATCC_17802_Vibrio_parahaemolyticus ATCC_19119_Listeria_ivanovii ATCC_25922_Escherichia_coli ATCC_33560_Campylobacter_jejuni ATCC_35221_Campylobacter_lari ATCC_35897_Listeria_welshimeri; do
    compare_assemblies.py --padding 50 --merge 100 --aligner edlib "$a".fasta.gz ../reference_assemblies/"$a".fasta.gz | paste - - - - | cut -f2 -d$'\t' | perl -pe 's/^ */>'$a'_/' | sed 's/: /\n/' | sed 's/ /_/' | sed 's/-//g' > "$a"_error_locations.fasta
done

# Align the error locations back to the references with bwa, aligning to all
# possible locations:
rm -f error_locations.sam
for a in ATCC_10708_Salmonella_enterica ATCC_14035_Vibrio_cholerae ATCC_17802_Vibrio_parahaemolyticus ATCC_19119_Listeria_ivanovii ATCC_25922_Escherichia_coli ATCC_33560_Campylobacter_jejuni ATCC_35221_Campylobacter_lari ATCC_35897_Listeria_welshimeri; do
    bwa index ../reference_assemblies/"$a".fasta.gz
    bwa mem -t 16 -a ../reference_assemblies/"$a".fasta.gz "$a"_error_locations.fasta >> error_locations.sam
    rm ../reference_assemblies/*.amb ../reference_assemblies/*.ann ../reference_assemblies/*.bwt ../reference_assemblies/*.pac ../reference_assemblies/*.sa
done

# Looking at the error_locations.sam file, I can see that only one of the
# errors was in a repeat (ATCC_35221_Campylobacter_lari position 491989), and
# that was right on the edge of a repeat, so probably not a problem for
# short-read alignment.

rm error_locations.sam *_error_locations.fasta
