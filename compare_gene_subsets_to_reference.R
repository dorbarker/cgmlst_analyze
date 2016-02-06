library(getopt)
library(parallel)

initial_options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script_name <- sub(file.arg.name, "", initial_options[grep(file.arg.name, initial_options)])
source_path <- paste(dirname(script_name),"partition_metrics.R", sep = "/")
source(source_path)

compare_partitions <- function(seed, n_genes, iterations, thresh, ref_thresholds) {
    # Calculate the Adjusted Wallace Coefficient and Cluster Cohesion
    # between the current subset and the reference

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
                         ref_thresholds = ref_thresholds,
                         thresh = i,
                         n_genes = n_genes
                        )


    }, mc.cores = cores)
#     l <- list()
#     for (i in 1:ncol(ref_thresholds)) {
#         l2 <- list()
#         for (x in seq_along(method_iterations)) {
#             l2[[x]] <- compare_partitions(seed = x,
#                                n_genes = n_genes,
#                                iterations = method_iterations,
#                                ref_thresholds = ref_thresholds,
#                                thresh = i
#             )
#         }
#         l[[i]] <- l2
#     }
    print(str(l))
    # kludge to convert to nice data.frame
    df <- do.call('rbind', do.call('rbind', l))
    df <- df[, colnames(df)[c(3, 1, 2, 4, 5, 6, 7)] ]

    df
}

load_data <- function(precomp_cluster_dir, cores) {

    extract_info <- function(x) {
        # Assumes CSV tables are named in the format of "NNN_clusters.csv"
        # where the Ns are integers representing the number of genes
        # used to calculate the clusters

        name <- basename(tools::file_path_sans_ext(x))
        pos <- c(gregexpr("\\d", name)[[1]])
        int_name <- as.integer(substr(name, min(pos), max(pos)))

        list(
            genes = int_name,
            clusters = read.csv(x, row.names = 1, check.names = FALSE)
            )

    }

    files <- list.files(precomp_cluster_dir, full.names = TRUE)
    reference_index <- grep("reference\\.csv", files)

    reference <- read.csv(files[reference_index],
                          row.names = 1,
                          check.names = FALSE)

    files <- files[-reference_index]

    list("submethods" = lapply(files, extract_info ),
         "reference" = reference)
}
options <- function() {

    spec <- matrix(c(
        "help",  "h", "0", "logical",   "Print this help and exit",
        "input", "i", "1", "character", "Path to directory of cluster tables",
        "out",   "o", "1", "character", "Output file",
        "cores", "c", "1", "integer",   "Number of CPU cores to use"
    ), byrow = TRUE, ncol = 5)

    opt <- getopt(spec)

    if (!is.null(opt$help)) {
        cat(getopt(spec, usage = TRUE))
        q(status = 1)
    }

    opt
}

main <- function() {

    opt <- options()
    data <- load_data(opt$input, opt$cores)

    l <- lapply(data$submethods, function(submethod) {

        compare_to_ref(n_genes = submethod$genes,
                       method_iterations = submethod$clusters,
                       ref_thresholds = data$reference,
                       cores = opt$cores)
    })
    print("2")
    out_df <- do.call('rbind', l)
    print("3")
    write.csv(out_df, file = opt$out, row.names = FALSE, quote = FALSE)
    print("4")
}

main()
