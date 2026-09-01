#!/bin/bash

# Activate the eggNOG environment
conda activate eggnog  

# Define paths
DATA_DIR="$(pwd)/eggnog-mapper/data"

# Check if the eggNOG database directory and core database file exist
if [[ ! -d "$DATA_DIR" ]]; then
    echo "Error: Database directory not found at $DATA_DIR"
    echo "Please create the directory or update DATA_DIR in this script."
    exit 1
fi

if [[ ! -f "$DATA_DIR/eggnog.db" ]]; then
    echo "Error: eggnog.db not found in $DATA_DIR"
    echo "Please ensure eggNOG-mapper databases are downloaded before running."
    exit 1
fi

echo "Database verified at $DATA_DIR"

# Create output directory
mkdir -p eggnog_annotations

# List of organisms
org_list=("I46" "I48" "YL2" "YL31" "YL32" "YL44" "YL45" "YL58" "KB1" "KB18")

# Loop over organisms
for org in "${org_list[@]}"; do
    echo "Processing $org..."

    cleaned_faa="fasta_protein/${org}_cleaned.faa"
    output_file="eggnog_annotations/${org}.emapper.annotations"

    # Verify input FASTA exists
    if [[ ! -f "$cleaned_faa" ]]; then
        echo "Warning: $cleaned_faa not found. Skipping $org..."
        continue
    fi

    # Run eggnog-mapper if output does not exist
    if [[ ! -f "$output_file" ]]; then
        emapper.py -i "$cleaned_faa" \
                   -o "eggnog_annotations/${org}" \
                   --cpu 16 \
                   --data_dir "$DATA_DIR"
        echo "eggnog-mapper finished for $org"
    else
        echo "Output already exists: $output_file, skipping..."
    fi

done

echo "All organisms processed."
