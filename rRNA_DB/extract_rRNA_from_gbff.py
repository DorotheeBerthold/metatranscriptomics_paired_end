#!/Users/dobertho/Desktop/Dorothee_homeoffice/R_analysis/metagenomics_tutorial/yes/bin/python
from Bio import SeqIO
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

with open(outfile, "w") as out:
    for record in SeqIO.parse(infile, "genbank"):
        for feature in record.features:
            if feature.type == "rRNA":
                product = feature.qualifiers.get("product", ["unknown_rRNA"])[0]
                product = product.replace(" ", "_").replace(",", "_").replace("(", "_").replace(")", "_")
                seq = feature.extract(record.seq)
                out.write(f">{record.id}_{product}\n{seq}\n")
