#!/bin/bash

# Exit on error, undefined variable, or pipeline failure
set -euo pipefail

# Script to create soft links to required files for the combine_predicted_novel_last_exons.py script
# Takes input and output directories as positional arguments

# Function to display usage information
display_help() {
    echo "Usage: $0 [input_dir] [output_dir]"
    echo
    echo "Create symbolic links to required files for combine_predicted_novel_last_exons.py across multiple PAPA identification runs for different datasets."
    echo
    echo "Arguments:"
    echo "  input_dir    Path to directory containing experiment subdirectories"
    echo "  output_dir   Path where links will be created (one subdirectory per experiment/subdirectory under input_dir)"
    echo
    echo "Required files (automatically linked/assumed to exist in experiment subdirectories under input_dir):"
    echo "  - tx_filtering/all_conditions.merged_last_exons.3p_end_filtered.gtf"
    echo "  - tx_filtering/novel_ref_combined.tx2le.tsv"
    echo "  - tx_filtering/novel_ref_combined.le2gene.tsv"
    echo "  - tx_filtering/novel_ref_combined.le2genename.tsv"
    echo "  - differential_apa/summarised_pas_quantification.ppau.tsv"
    echo
    echo "Example: $0 /path/to/experiments /path/to/links"
    exit 1
}

# Check if required arguments are provided
if [ $# -ne 2 ]; then
    echo "Error: Exactly two arguments required."
    display_help
fi

# Check if help is requested
if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    display_help
fi


# Set input and output directories
input_dir="$1"
output_dir="$2"

# Check if input directory exists
if [ ! -d "$input_dir" ]; then
    echo "Error: Input directory '$input_dir' does not exist."
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each subdirectory in the input directory
for dir in "$input_dir"/*/; do 
    # Get basename of directory (remove trailing slash)
    base_dir=$(basename "${dir%/}")
    
    echo "Processing directory: $base_dir"
    
    # Create subdirectory in output directory
    mkdir -p "$output_dir/$base_dir"
    
    # Create symbolic links to required files
    ln -s "$dir/tx_filtering/all_conditions.merged_last_exons.3p_end_filtered.gtf" "$output_dir/$base_dir/" 
    ln -s "$dir/tx_filtering/novel_ref_combined.tx2le.tsv" "$output_dir/$base_dir/" 
    ln -s "$dir/tx_filtering/novel_ref_combined.le2gene.tsv" "$output_dir/$base_dir/" 
    ln -s "$dir/tx_filtering/novel_ref_combined.le2genename.tsv" "$output_dir/$base_dir/" 
    ln -s "$dir/differential_apa/summarised_pas_quantification.ppau.tsv" "$output_dir/$base_dir/"
    
    echo "  Created links for $base_dir"
done

echo "Done! Links created in $output_dir"
