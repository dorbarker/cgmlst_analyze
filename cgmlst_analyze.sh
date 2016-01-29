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
            echo "--cores            Number of CPU cores to use for parallel tasks (default = all available)."
            echo "--gene-select      Comma-seperated list of numbers of genes to subset."
            echo "--replicates       Number of replicates for gene subsetting (default = 1000)."
            echo "--prediction-width Comma-seperated list of number of genes to consider for allele prediction."
            exit 0
            ;;
        
        --work-dir)
            shift
            if test $# -gt 0; then
                export WORKDIR=$(echo $(fullpath "$1"))/
            else
                echo "You need to give me a work directory."
                exit 1
            fi
            shift
            ;;

        --cores)
            shift
            if test $# -gt 0; then
                export CORES=$1
            else
                # Will use all available cores if not specified
                export CORES=$( grep -ci processor /proc/cpuinfo )
            fi
            shift
            ;;

        --gene-select)
            shift
            if test $# -gt 0; then
                export GENESELECT=$1
            else
                echo "You must provide a comma-seperated list of gene counts to subset"
                exit 1
            fi
            shift
            ;;

        --replicates)
            shift
            if test $# -gt 0; then
                export REPLICATES=$1
            else
                export REPLICATES=1000
            fi
            shift
            ;;
        --prediction-widths)
            shift
            if test $# -gt 0; then
               export PREDWIDTH=$1
            else
                echo "You need to specify allele prediction widths."
            fi
            shift
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


### Generate gene subsets ###
printf "\nGenerating gene subsets.\n"
Rscript ${SCRIPT_DIR}/pre_compute_clusters.R --input core_calls.csv \
                                             --n_select $GENESELECT \
                                             --replicates $REPLICATES\
                                             --cores $CORES \
                                             --outdir ${WORKDIR}/pre_computed_clusters/

### Compare gene subsets to reference clusters ###
printf "\nComparing gene subsets to the reference set.\n"
Rscript ${SCRIPT_DIR}/compare_gene_subsets_to_reference.R --input ${WORKDIR}/pre_computer_clusters \
                                                          --out   ${WORKDIR}/gene_subset_metrics.csv \
                                                          --cores $CORES

# Optimized subregions
    # Map contigs
    # Shannon entropy
        # k-mer
        # base
        # gene
    # Genetic Algorithm 
    # Naive subregions

# Allele imputation
    
### Naive association ###
printf "\nFinding allele associations - naive.\n"
mkdir ${WORKDIR}/triplets
Rscript ${SCRIPT_DIR}/triplet_linkage.R --input ${WORKDIR}/core_calls.csv \
                                        --out  ${WORKDIR}/triplets/ \
                                        --cores $CORES
### Machine learning ###
printf "\nFinding allele associations - random forest\n"
Rscript ${SCRIPT_DIR}/allele_associations_ml.R  --input ${WORKDIR}/core_calls.csv \
                                                --out ${WORKDIR}/ \
                                                --min 5 \
                                                --width $PREDWIDTH \
                                                --train 0.1 \
                                                --trees $REPLICATES \
                                                --cores $CORES


