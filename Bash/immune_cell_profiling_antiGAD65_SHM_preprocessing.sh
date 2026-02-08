#!/bin/bash
################################################################################
# BCR Repertoire Processing Pipeline using Immcantation Suite
################################################################################
#
# Description: Process 10X Genomics single-cell BCR sequencing data using
#              the Immcantation Suite to identify clones and reconstruct
#              germline sequences in AIRR Community standardized format.
#
# Author: Sumanta Barman
# Date: February 2026
# Version: 1.0
#
# Reference: Immcantation Suite v4.5.0
#            https://immcantation.readthedocs.io/
#
# Requirements:
#   - Docker or Singularity
#   - 10X Genomics VDJ output files:
#       * filtered_contig.fasta
#       * filtered_contig_annotations.csv
#
# Output: AIRR-formatted database with clonal assignments
#
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

################################################################################
# CONFIGURATION
################################################################################

# Docker/Singularity settings
CONTAINER_IMAGE="immcantation/suite:4.5.0"
ORGANISM="human"
LOCI="ig"

# Input files (modify these paths as needed)
INPUT_FASTA="filtered_contig.fasta"
INPUT_ANNOTATIONS="filtered_contig_annotations.csv"
OUTPUT_DIR="results"
OUTPUT_NAME="BCR_data_sequences"

# Clonal clustering parameters
DISTANCE_THRESHOLD="0.16"  # Normalized Hamming distance threshold

# Germline reference paths (inside container)
IGBLAST_DB="/usr/local/share/igblast"
GERMLINE_DIR="/usr/local/share/germlines/imgt/human/vdj"

################################################################################
# FUNCTIONS
################################################################################

# Print colored messages
print_header() {
    echo ""
    echo "================================================================================"
    echo "$1"
    echo "================================================================================"
    echo ""
}

print_step() {
    echo ""
    echo ">>> Step $1: $2"
    echo ">>> $(date '+%Y-%m-%d %H:%M:%S')"
    echo ""
}

print_success() {
    echo "✓ $1"
}

# Check if required files exist
check_files() {
    print_header "Checking Input Files"

    if [ ! -f "$INPUT_FASTA" ]; then
        echo "ERROR: Input FASTA file not found: $INPUT_FASTA"
        exit 1
    fi
    print_success "Found: $INPUT_FASTA"

    if [ ! -f "$INPUT_ANNOTATIONS" ]; then
        echo "ERROR: Input annotations file not found: $INPUT_ANNOTATIONS"
        exit 1
    fi
    print_success "Found: $INPUT_ANNOTATIONS"

    echo ""
}

################################################################################
# MAIN PIPELINE
################################################################################

print_header "BCR Repertoire Processing Pipeline - Immcantation Suite v4.5.0"

# Check input files
check_files

# Pull Docker image
print_step "0" "Pulling Docker Image"
docker pull $CONTAINER_IMAGE
print_success "Docker image pulled successfully"

# Optional: Start interactive container (for manual exploration)
# Uncomment the following line to run interactively:
# docker run -it $CONTAINER_IMAGE bash

################################################################################
# STEP 1: ASSIGN V(D)J GENES AND CONVERT TO AIRR FORMAT
################################################################################

print_step "1" "Assigning V(D)J Genes using IgBLAST"

docker run --rm \
    -v "$(pwd):/data" \
    -w /data \
    $CONTAINER_IMAGE \
    AssignGenes.py igblast \
        -s "$INPUT_FASTA" \
        -b "$IGBLAST_DB" \
        --organism "$ORGANISM" \
        --loci "$LOCI" \
        --format blast \
        --outdir "$OUTPUT_DIR" \
        --outname "$OUTPUT_NAME"

print_success "V(D)J gene assignment complete"
print_success "Output: ${OUTPUT_DIR}/${OUTPUT_NAME}_igblast.fmt7"

################################################################################
# STEP 2: CREATE AIRR-FORMATTED DATABASE
################################################################################

print_step "2" "Creating AIRR-Formatted Database"

docker run --rm \
    -v "$(pwd):/data" \
    -w /data \
    $CONTAINER_IMAGE \
    MakeDb.py igblast \
        -i "${OUTPUT_DIR}/${OUTPUT_NAME}_igblast.fmt7" \
        -s "$INPUT_FASTA" \
        -r "$GERMLINE_DIR" \
        --10x "$INPUT_ANNOTATIONS" \
        --extended

print_success "AIRR database created"
print_success "Output: ${OUTPUT_DIR}/${OUTPUT_NAME}_igblast_db-pass.tsv"

################################################################################
# STEP 3: SEPARATE HEAVY AND LIGHT CHAINS
################################################################################

print_step "3a" "Parsing Heavy Chain Sequences (IGH)"

docker run --rm \
    -v "$(pwd):/data" \
    -w /data \
    $CONTAINER_IMAGE \
    ParseDb.py select \
        -d "${OUTPUT_DIR}/${OUTPUT_NAME}_igblast_db-pass.tsv" \
        -f locus \
        -u "IGH" \
        --logic all \
        --regex \
        --outname BCR_heavy

print_success "Heavy chain sequences extracted"
print_success "Output: ${OUTPUT_DIR}/BCR_heavy_parse-select.tsv"

print_step "3b" "Parsing Light Chain Sequences (IGK/IGL)"

docker run --rm \
    -v "$(pwd):/data" \
    -w /data \
    $CONTAINER_IMAGE \
    ParseDb.py select \
        -d "${OUTPUT_DIR}/${OUTPUT_NAME}_igblast_db-pass.tsv" \
        -f locus \
        -u "IG[LK]" \
        --logic all \
        --regex \
        --outname BCR_light

print_success "Light chain sequences extracted"
print_success "Output: ${OUTPUT_DIR}/BCR_light_parse-select.tsv"

################################################################################
# STEP 4: DEFINE CLONES FROM HEAVY CHAIN
################################################################################

print_step "4" "Clustering Heavy Chain Sequences into Clones"

docker run --rm \
    -v "$(pwd):/data" \
    -w /data \
    $CONTAINER_IMAGE \
    DefineClones.py \
        -d "${OUTPUT_DIR}/BCR_heavy_parse-select.tsv" \
        --act set \
        --model ham \
        --norm len \
        --dist "$DISTANCE_THRESHOLD"

print_success "Clonal clustering complete"
print_success "Output: ${OUTPUT_DIR}/BCR_heavy_parse-select_clone-pass.tsv"
echo "    Distance threshold: $DISTANCE_THRESHOLD (normalized Hamming distance)"

################################################################################
# STEP 5: ASSIGN LIGHT CHAINS TO HEAVY CHAIN CLONES
################################################################################

print_step "5" "Assigning Light Chains to Heavy Chain Clones"

docker run --rm \
    -v "$(pwd):/data" \
    -w /data \
    $CONTAINER_IMAGE \
    light_cluster.py \
        -d "${OUTPUT_DIR}/BCR_heavy_parse-select_clone-pass.tsv" \
        -e "${OUTPUT_DIR}/BCR_light_parse-select.tsv" \
        -o "${OUTPUT_DIR}/BCR_heavy_light_clone_parse.tsv"

print_success "Light chains assigned to clones"
print_success "Output: ${OUTPUT_DIR}/BCR_heavy_light_clone_parse.tsv"

################################################################################
# STEP 6: RECONSTRUCT GERMLINE SEQUENCES
################################################################################

print_step "6" "Reconstructing Germline Sequences"

docker run --rm \
    -v "$(pwd):/data" \
    -w /data \
    $CONTAINER_IMAGE \
    CreateGermlines.py \
        -d "${OUTPUT_DIR}/BCR_heavy_parse-select_clone-pass.tsv" \
        -g dmask \
        --cloned \
        -r "${GERMLINE_DIR}/IGHV.fasta" \
           "${GERMLINE_DIR}/IGHD.fasta" \
           "${GERMLINE_DIR}/IGHJ.fasta"

print_success "Germline reconstruction complete"
print_success "Output: ${OUTPUT_DIR}/BCR_heavy_parse-select_clone-pass_germ-pass.tsv"

################################################################################
# PIPELINE COMPLETION
################################################################################

print_header "Pipeline Complete!"

echo "Summary of Output Files:"
echo "------------------------"
echo "1. IgBLAST results:     ${OUTPUT_DIR}/${OUTPUT_NAME}_igblast.fmt7"
echo "2. AIRR database:       ${OUTPUT_DIR}/${OUTPUT_NAME}_igblast_db-pass.tsv"
echo "3. Heavy chain parsed:  ${OUTPUT_DIR}/BCR_heavy_parse-select.tsv"
echo "4. Light chain parsed:  ${OUTPUT_DIR}/BCR_light_parse-select.tsv"
echo "5. Clonal assignment:   ${OUTPUT_DIR}/BCR_heavy_parse-select_clone-pass.tsv"
echo "6. Heavy+Light clones:  ${OUTPUT_DIR}/BCR_heavy_light_clone_parse.tsv"
echo "7. Germline sequences:  ${OUTPUT_DIR}/BCR_heavy_parse-select_clone-pass_germ-pass.tsv"
echo ""
echo "Pipeline completed successfully at $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

################################################################################
# END OF PIPELINE
################################################################################
