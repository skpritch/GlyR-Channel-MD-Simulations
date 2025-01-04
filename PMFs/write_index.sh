#!/bin/bash

# Input and output directories
input_dir="david_CONFIGS_CLA"
output_dir="DAVID_index_CLA"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop over numbers from -45 to 45
for num in $(seq -45 45); do
    # Construct input and output file names
    gro_file="$input_dir/conf${num}.gro"
    ndx_file="$output_dir/index${num}.ndx"

    # Check if the .gro file exists
    if [ -f "$gro_file" ]; then
        echo "Processing $gro_file..."

        # Run make_ndx
        gmx make_ndx -f "$gro_file" -o "$ndx_file" << EOF
a 28728
name 17 ion
r POPC
name 18 MEMB
r TIP3 | r POT | r CLA
name 19 SOLV
r ALA | r ARG | r ASN | r ASP | r CYS | r GLN | r GLU | r GLY | r HSD | r ILE | r LEU | r LYS | r MET | r PHE | r PRO | r SER | r THR | r TRP | r TYR | r VAL
name 20 SOLU
"SOLU" | "MEMB"
name 21 SOLU_MEMB
q
EOF

    else
        echo "File $gro_file not found. Skipping..."
    fi
done

echo "Index file generation completed. Index files are in the directory '$output_dir'."

