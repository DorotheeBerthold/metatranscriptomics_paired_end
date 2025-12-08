#!/Users/dobertho/Desktop/Dorothee_homeoffice/R_analysis/metagenomics_tutorial/yes/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  2 14:20:37 2025

@author: dobertho
"""


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation
import os

gbff_dir = "OMM_gbff"
pseudogenome_fasta = "OMM12_pseudogenome.fasta"
pseudogenome_gff = "OMM12_pseudogenome.gff3"

with open(pseudogenome_fasta, "w") as fasta_out, open(pseudogenome_gff, "w") as gff_out:
    gff_out.write("##gff-version 3\n")
    
    for gbff_file in sorted(os.listdir(gbff_dir)):
        if not gbff_file.endswith(".gbff"):
            continue

        strain = gbff_file.replace(".gbff", "")
        record = SeqIO.read(os.path.join(gbff_dir, gbff_file), "genbank")
	
        # Write strain sequence
        record.id = strain
        record.description = ""
        record.name = strain
        SeqIO.write(record, fasta_out, "fasta")

        # Write features
        for feature in record.features:
            if feature.type == "source":
                continue
            loc = feature.location

    # Handle join() / circular wrap-around
            if isinstance(loc, CompoundLocation):
                start = int(loc.parts[0].start) + 1
                end   = int(loc.parts[-1].end)
            else:
                start = int(loc.start) + 1
                end   = int(loc.end)

# IMPORTANT FIX:
            strand_val = loc.strand

    # Correct way to get strand
            if strand_val is None:
                strand = "."
            elif strand_val == 1:
                strand = "+"
            elif strand_val == -1:
                strand = "-"
            else:
                strand = "."

            feature_type = feature.type

            attributes = []
            for key in ["gene", "product", "locus_tag", "note"]:
                if key in feature.qualifiers:
                    attributes.append(f"{key}={feature.qualifiers[key][0]}")
            attr_str = ";".join(attributes) if attributes else "."

            gff_out.write(f"{strain}\tOMM12\t{feature_type}\t{start}\t{end}\t.\t{strand}\t.\t{attr_str}\n")
