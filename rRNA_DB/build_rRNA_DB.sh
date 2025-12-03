# Required packages: infernal & mafft
#First activate extract rRNA_from_gbff.py
chmod +x extract_rRNA_from_gbff.py


# Run extract all rRNA from files
chmod +x run_extract_all.sh
./run_extract_all.sh

# Split the files according to rRNA type
grep -A1 -E "16S|SSU" OMM12_rRNAs.fasta > OMM12_16S.fasta
grep -A1 -E "23S|LSU" OMM12_rRNAs.fasta > OMM12_23S.fasta
grep -A1 -E "5S" OMM12_rRNAs.fasta  > OMM12_5S.fasta

# Perform multiple sequence alignement using MAFFT
mafft --auto OMM12_16S.fasta > OMM12_16S.aln.fasta
mafft --auto OMM12_23S.fasta > OMM12_23S.aln.fasta
mafft --auto OMM12_5S.fasta  > OMM12_5S.aln.fasta


# Convert to Stockholm format
esl-reformat stockholm OMM12_16S.aln.fasta > OMM12_16S.sto
esl-reformat stockholm OMM12_23S.aln.fasta > OMM12_23S.sto
esl-reformat stockholm OMM12_5S.aln.fasta  > OMM12_5S.sto


# Build covariance models using Infernal ignoring secondary structures

cmbuild --noss OMM12_16S.cm OMM12_16S.sto
cmbuild --noss OMM12_23S.cm OMM12_23S.sto
cmbuild --noss OMM12_5S.cm  OMM12_5S.sto


# Combine all covariance models into a single database
cat OMM12_16S.cm OMM12_23S.cm OMM12_5S.cm > OMM12_rRNA.cm
# Press CM database
cmpress OMM12_rRNA.cm
