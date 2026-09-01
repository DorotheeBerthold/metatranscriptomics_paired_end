#!/bin/bash

# Define your paths here to keep the script clean
DIR_INPUT="cleaned_faa files"
DIR_OUTPUT="cleaned_faa files/cayman"
DB_PATH="v3/" 
CUTOFFS="v3/cutoffs.csv" 

mkdir -p "$DIR_OUTPUT"

# --- 2. PROCESSING LOOP ---
# Loop through every .faa file in the input directory
for INPUT_PATH in "$DIR_INPUT"/*.faa; do
    
    # Handle case where no .faa files are found
    [[ -e "$INPUT_PATH" ]] || continue

    # Extract just the filename (e.g., "sample.faa")
    FILENAME=$(basename "$INPUT_PATH")
    
    # Extract the name without extension (e.g., "sample")
    BASENAME="${FILENAME%.*}"

    # Define the output file path
    OUTPUT_FILE="$DIR_OUTPUT/${BASENAME}_cayman.csv"

    # Check if output file already exists
    if [[ -f "$OUTPUT_FILE" ]]; then
        echo "------------------------------------------------"
        echo "Output already exists for $FILENAME, skipping..."
        continue
    fi

    echo "------------------------------------------------"
    echo "Processing: $FILENAME"

    # --- 3. RUN CAYMAN ---
    cayman annotate_proteome \
        "$DB_PATH" \
        "$INPUT_PATH" \
        -o "$OUTPUT_FILE" \
        --cutoffs "$CUTOFFS" \
        --threads 4

    echo "Finished: $BASENAME"
done

echo "------------------------------------------------"
echo "All files processed successfully."
