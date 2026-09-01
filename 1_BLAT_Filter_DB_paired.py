#!/usr/bin/env python

import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator

if len(sys.argv) != 9:
    print("Usage: ./1_BLAT_Filter_DB_paired.py <R1_in.fq> <R2_in.fq> <R1.blat> <R2.blat> <R1_clean.fq> <R2_clean.fq> <R1_contam.fq> <R2_contam.fq>")
    sys.exit(1)

r1_in_file = sys.argv[1]
r2_in_file = sys.argv[2]
blat_r1_file = sys.argv[3]
blat_r2_file = sys.argv[4]
r1_clean_out = sys.argv[5]
r2_clean_out = sys.argv[6]
r1_contam_out = sys.argv[7]
r2_contam_out = sys.argv[8]

# 1. Gather all contaminated read IDs from BOTH BLAT files
# If a read hits the host database in either R1 or R2, we discard the whole pair.
contaminated = set()
for blat_file in [blat_r1_file, blat_r2_file]:
    with open(blat_file) as f:
        for line in f:
            if len(line.strip()) > 0:
                query_id = line.split("\t")[0]
                contaminated.add(query_id)

output_count = 0
contaminant_count = 0

# 2. Process both FASTQ files simultaneously to keep pairs perfectly synced
with open(r1_in_file) as f1_in, open(r2_in_file) as f2_in, \
     open(r1_clean_out, "w") as f1_clean, open(r2_clean_out, "w") as f2_clean, \
     open(r1_contam_out, "w") as f1_contam, open(r2_contam_out, "w") as f2_contam:

    # zip() lets us step through R1 and R2 side-by-side
    for (id1, seq1, qual1), (id2, seq2, qual2) in zip(FastqGeneralIterator(f1_in), FastqGeneralIterator(f2_in)):
        
        # The ID in BLAT is usually just the first part of the header (before the space)
        base_id = id1.split(" ")[0]

        if base_id in contaminated:
            # Write to contaminated files
            f1_contam.write(f"@{id1}\n{seq1}\n+\n{qual1}\n")
            f2_contam.write(f"@{id2}\n{seq2}\n+\n{qual2}\n")
            contaminant_count += 1
        else:
            # Write to clean files
            f1_clean.write(f"@{id1}\n{seq1}\n+\n{qual1}\n")
            f2_clean.write(f"@{id2}\n{seq2}\n+\n{qual2}\n")
            output_count += 1

print(f"{contaminant_count} read pairs were aligned to the contaminant database")
print(f"{output_count} read pairs were completely clean")