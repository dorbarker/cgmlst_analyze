library(getopt)
library(parallel)

source("./partition_metrics.R")

compare_partitions <- function(seed, n_genes, iterations, thresh, ref_thresholds) {

    iteration <- iterations[,seed]
    ref_thresh <- ref_thresholds[,thresh]


    awc <- adj_wallace(iteration, ref_thresh)

    list(
    seed = seed,
    reference_threshold = as.numeric(colnames(ref_thresholds)[thresh]),
    genes = n_genes,
    AWC_subset_vs_ref = awc$Adjusted_Wallace_A_vs_B,
    AWC_ref_vs_subset = awc$Adjusted_Wallace_B_vs_A,

    cluster_cohesion_subset_vs_ref = taboada_cohesion(iteration, ref_thresh)$weighted_average,
    cluster_cohesion_ref_vs_subset = taboada_cohesion(ref_thresh, iteration)$weighted_average
    )
}

compare_to_ref <- function(n_genes, method_iterations, ref_thresholds, cores) {


    l <- mclapply(1:ncol(ref_thresholds), function(i) {

        lapply(seq_along(method_iterations),
                         compare_partitions,
                         iterations = method_iterations,
                         ref_thresh = ref_thresholds,
                         thresh = i,
                         n_genes = n_genes
                        )


    }, mc.cores = cores)

    # kludge to convert to nice data.frame
    df <- do.call('rbind', do.call('rbind', l))
    df <- df[, colnames(df)[c(3, 1, 2, 4, 5, 6, 7)] ]
}

load_data <- function(precomp_cluster_dir, cores) {

    extract_info <- function(x) {

        name <- basename(tools::file_path_sans_ext(x))
        pas <- c(gregexpr("\\d", name))
        int_name <- as.integer(substr(name, min(pos), max(pos)))

        list(
            genes = int_name,
            clusters = read.csv(x, row.names = 1, check.names = FALSE)
            )

    }

    files <- list.files(precomp_cluster_dir, full.names = TRUE)
    reference_index <- grep(files, "reference\\.csv")

    reference <- read.csv(files[reference_index],
                          row.names = 1,
                          check.names = FALSE)

    files <- files[-reference_index]

    list("submethods" = mclapply(files, extract_info, mc.cores = cores),
         "reference" = reference)

}

main <- function() {

    opt <- options()



}
