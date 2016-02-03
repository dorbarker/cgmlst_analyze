library(parallel)
library(ape)
library(getopt)

get_ordered_clusters <- function(ref_calls, selection) {

    # Find 100% clusters, and return them sorted by strain name

    d <- dist.gene(ref_calls[, selection])
    hc <- hclust(d, "complete")

    clusts <- cutree(hc, h = 0)

    clusts[order(names(clusts))]
}

ref_clusters <- function(ref_calls, outdir) {

    d <- dist.gene(ref_calls)
    hc <- hclust(d, "complete")

    # Cuts at all possible numbers of gene differences
    clusts <- cutree(hc, h = unique(c(0, 0:ncol(ref_calls))))

    reference_df <- data.frame(clusts[order(names(rowclusts)),],
                               stringsAsFactors = FALSE,
                               check.names = FALSE)

    write.csv(reference_df, file = paste0(outdir, "/", "reference.csv"),
              row.names = TRUE, quote = FALSE)
}

bootstrap_factory <- function(n_select, ref_calls) {

    n_genes <- ncol(ref_calls)

    f <- function(seed) {
        # Closure to be used in compute_clusters()'s  mclapply for simplicity

        # sets the seed from the input argument of mclapply
        set.seed(seed)

        # subsamples genes from ref_calls
        selection <- sample(1:n_genes, size = n_select, replace = FALSE)

        # clusters from selected genes
        clusters <- get_ordered_clusters(ref_calls, selection)

    }

    f

}

compute_clusters <- function(n_select, ref_calls, cores, replicates) {

    f <- bootstrap_factory(n_select, ref_calls)

    # Setting seed allows both reproducibility
    # and nesting of gene selections
    seeds <- 1:replicates
    clusters <- mclapply(seeds, f, mc.cores = cores, mc.set.seed = FALSE)

    cluster_df <- as.data.frame(clusters, stringsAsFactors = FALSE,
                                check.names = FALSE)
    colnames(cluster_df) <- seeds

    cluster_df
}

write_clusters <- function(n_select, cluster_df, outdir) {

    # Create the output directory if it does not exist
    if (!dir.exists(outdir)) {
        dir.create(outdir)
    }

    outname <- paste0(outdir, "/", n_select, "_clusters.csv")

    write.csv(cluster_df, file = outname, row.names = TRUE, quote = FALSE)

    cat(paste("\n", "Clusters written to:", outname, "\n"))
}

options <- function() {

    n_select_help <- paste("Number of genes to select;",
                           "multiple values may be given",
                           "as comma-seperated list")
    # commandline arguments
    spec <- matrix(c(
        "help", "h", "0", "logical", "Print this help and exit",
        "input", "i", 1, "character", "Input allele calls",
        "n_select", "n", 1, "character", n_select_help,
        "replicates", "r", 1, "integer", "Number of replicates",
        "cores", "c", 1, "integer", "Number of CPU cores to use",
        "outdir", "o", 1, "character", "Output path for clusters"
    ), byrow = TRUE, ncol = 5)

    opt <- getopt(spec)

    if (!is.null(opt$help)) {
        cat(getopt(spec, usage = TRUE))
        q(status = 1)
    }
    opt$n_select <- as.integer(strsplit(opt$n_select, ","))

    opt
}

main <- function() {

    opt <- options()

    # Load reference calls from which subsets will be chosen
    ref_calls <- read.csv(opt$input, stringsAsFactors = FALSE, row.names = 1)

    # Computer ref clusters over all thresholds and write out
    ref_clusters(ref_calls, opt$outdir)

    for (n_gene in opt$n_select) {
    cluster_df <- compute_clusters(n_select = n_gene,
                                   ref_calls = ref_calls,
                                   cores = opt$cores,
                                   replicates = opt$reps)

    write_clusters(opt$n_select, cluster_df, opt$outdir)
    }
}

main()
