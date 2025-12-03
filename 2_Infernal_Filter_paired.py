#!/usr/bin/env python

import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

r1_file = sys.argv[1]
r2_file = sys.argv[2]
infernal_out = sys.argv[3]
mRNA_r1_file = sys.argv[4]
mRNA_r2_file = sys.argv[5]
rRNA_r1_file = sys.argv[6]
rRNA_r2_file = sys.argv[7]

# Load Infernal IDs
Infernal_rRNA_IDs = set()
with open(infernal_out) as f:
    for line in f:
        if not line.startswith("#"):
            Infernal_rRNA_IDs.add(line.split()[0])

# Open files
r1_iter = SeqIO.parse(r1_file, "fastq")
r2_iter = SeqIO.parse(r2_file, "fastq")

mRNA_r1 = open(mRNA_r1_file, "w")
mRNA_r2 = open(mRNA_r2_file, "w")
rRNA_r1 = open(rRNA_r1_file, "w")
rRNA_r2 = open(rRNA_r2_file, "w")

mRNA_r1_count = mRNA_r2_count = 0
rRNA_r1_count = rRNA_r2_count = 0
discarded = 0

for read1, read2 in zip(r1_iter, r2_iter):

    id1 = read1.id.split()[0]
    id2 = read2.id.split()[0]

    if id1 != id2:
        print("Mismatch detected")
        print("R1:", read1.id)
        print("R2:", read2.id)
        sys.exit(1)

    # Read IDs must match
    assert read1.id.split()[0] == read2.id.split()[0]

    in_rRNA1 = read1.id in Infernal_rRNA_IDs
    in_rRNA2 = read2.id in Infernal_rRNA_IDs
    
    if read1.id.split()[0] != read2.id.split()[0]:
        print("Mismatch:")
        print("R1:", read1.id)
        print("R2:", read2.id)
        sys.exit(1)

    if in_rRNA1 and in_rRNA2:
        SeqIO.write(read1, rRNA_r1, "fastq")
        SeqIO.write(read2, rRNA_r2, "fastq")
        rRNA_r1_count += 1
    elif not in_rRNA1 and not in_rRNA2:
        SeqIO.write(read1, mRNA_r1, "fastq")
        SeqIO.write(read2, mRNA_r2, "fastq")
        mRNA_r1_count += 1
    else:
        # One mate is rRNA, the other is mRNA → discard pair
        discarded += 1

print(f"mRNA pairs kept: {mRNA_r1_count}")
print(f"rRNA pairs kept: {rRNA_r1_count}")
print(f"discordant pairs discarded: {discarded}")
