#!/bin/bash

# Activate the eggNOG environment
conda activate eggnog  

# List of organisms
org_list=("I46" "I48" "YL2" "YL31" "YL32" "YL44" "YL45" "YL58")

# Loop over organisms
for org in "${org_list[@]}"; do
    echo "Processing $org..."

    # Clean FASTA if it doesn't exist
    cleaned_faa="fasta_protein/${org}_cleaned.faa"
    if [[ ! -f $cleaned_faa ]]; then
        sed 's/ //g' "fasta_protein/${org}_proteins.faa" > "$cleaned_faa"
        echo "Cleaned FASTA created: $cleaned_faa"
    else
        echo "Cleaned FASTA already exists: $cleaned_faa"
    fi

    # Run eggnog-mapper if output does not exist
    output_file="eggnog_annotations/${org}.emapper.annotations"
    if [[ ! -f $output_file ]]; then
        emapper.py -i "$cleaned_faa" -o "eggnog_annotations/${org}" --cpu 16
        echo "eggnog-mapper finished for $org"
    else
        echo "Output already exists: $output_file, skipping..."
    fi

done

echo "All organisms processed."
