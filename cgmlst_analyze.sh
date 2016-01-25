#!/bin/bash


if test $# -eq 0; then
    echo "No arguments provided!"
    exit 1
fi

SRC_DIR=$(cd "$(dirname "$0")/.."; pwd)
SCRIPT_DIR=$(cd "$(dirname "$0")"; pwd)
T="$(date +%s)"
function fullpath {

    echo "$( cd "$(dirname "$1")"; pwd)/$(basename "$1")"
}


while test $# -gt 0; do

    case "$1" in
        
        -h|--help)
            echo "Required arguments:"
            echo "--work-dir         Working directory for this script."
            echo "--reference        Fasta filename (not path) within --genomes path to be used as reference."
            echo "--genomes          Path to directory containing genomes as FASTAS."
            exit 0
            ;;

    esac
done

##### Analyses #####

# Subsampling of genes
    # AWC
    # Taboada Cohesion
    # SID
    # N clusters
    # Export computed clusters

# Optimized subregions
    # Map contigs
    # Shannon entropy
        # k-mer
        # base
        # gene
    # Antfarm optimization?
    # Naive subregions

# Allele imputation
    # Machine learning
    # Naive association
