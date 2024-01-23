for a in ATCC_10708_Salmonella_enterica ATCC_14035_Vibrio_cholerae ATCC_17802_Vibrio_parahaemolyticus ATCC_19119_Listeria_ivanovii ATCC_25922_Escherichia_coli ATCC_33560_Campylobacter_jejuni ATCC_35221_Campylobacter_lari ATCC_35897_Listeria_welshimeri; do
    compare_assemblies.py --aligner edlib "$a".fasta.gz ../reference_assemblies/"$a".fasta.gz > "$a".errors
done
