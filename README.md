# Algorithms

## Gene Subset Selection & Pre-Computation of Clusters

`pre_compute_clusters.R`

  1. Load the reference calls (the full set of identified cgMLST genes)
  2. Calculate reference clusters over all possible thresholds
  3. For each *N* genes and for each replicate:
    
    1. Replicate number is used as the seed for random selection
    2. Select from the reference calls *N* random genes
    3. Calculate complete-linkage clusters from selected genes with `hclust()`

  4. Write seeds and associated clusters to CSV

## Gene Subset Comparison

`compare_gene_subsets_to_reference.R`

  1. Load gene subsets from pre-computed cluster directory
  2. Load reference clusters for each threshold
  3. For each combination of gene subset and reference threshold:
    
    1. Calculate bi-directional Adjusted Wallace Coefficients
    2. Calculate Cluster Cohesion metric bi-directionally

  4. Write results to CSV

## Partial Match Analysis

`partial_match_analysis.py`

  1. Gather all observed alleles for all genes in MIST's output JSONs
  2. For each gene and for each replicate:

    1. Randomly select an allele
    2. In that allele, randomly select a pivot and a side
    3. Take the substring on that side of the pivot
    4. If the substring versus NCTC11168 genome has `evalue <= 10.0`, continue, else repick a pivot and side
    5. Return the reciprocal number of alleles to which the substring identically matches
  
  3. Write values to a CSV


