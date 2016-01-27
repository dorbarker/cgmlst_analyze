library(parallel)
library(getopt)

process_triples <- function(start_index, func, df) {

    func(df[,start_index:(start_index + 2)])

}

get_triples_by_flank <- function(df) {

    triples <- list()

    for (i in 1:nrow(df)) {
        line <- df[i,]

        flank <- paste(line[1], line[3], sep = "_")
        centre <- as.character(line[2])

        if (flank %in% names(triples)) {

            if (centre %in% names(triples[[flank]])) {
                triples[[flank]][[centre]] <- triples[[flank]][[centre]] + 1
            }
            else {
                triples[[flank]][[centre]] <- 1
            }

        }
        else {

            triples[[flank]] <- list()
            triples[[flank]][[centre]] <- 1
        }
    }

    triples
}

get_triples_by_centre <- function(df) {

    triples <- list()
    for (i in 1:nrow(df)) {

        line <- df[i,]
        flank <- paste(line[1], line[3], sep = "_")
        centre <- as.character(line[2])

        if (centre %in% names(triples)) {
            if (flank %in% names(triples[[centre]])) {
                triples[[centre]][[flank]] <- triples[[centre]][[flank]] + 1
            }
            else {
                triples[[centre]][[flank]] <- 1
            }
        }
        else {
            triples[[centre]] <- list()
            triples[[centre]][[flank]] <- 1
        }
    }

    triples
}

extract_stats <- function(triple_list) {

    largest_proportion <- function(name) {

        ul <- unlist(triple_list[[name]])
        max(ul) / sum(ul)
    }

    weighted_proportion <- function(select_names) {

        largest <- sapply(select_names, function(x) max(unlist(triple_list[[x]])))

        total <- sapply(select_names, function(x) sum(unlist(triple_list[[x]])))

        if (length(select_names) > 0) {

            sum(largest) / sum(total)
        }
        else {
            1
        }
    }

    multi_partner <- names(triple_list)[which(sapply(triple_list, length) > 1)]

    multipartners_n_partner <- sapply(multi_partner, function(x) length(triple_list[[x]]))

    non_single_names <- names(triple_list)[which(sapply(triple_list, function(x) sum(unlist(x))) > 1)]
    non_single <- sapply(non_single_names, function(x) sum(unlist(triple_list[[x]])))
    names(non_single) <- non_single_names

    list(
        "non_singletons" = non_single_names,
        "n_unique"  = length(unique(names(triple_list))),
        "genes_w_multiple_partners" = multi_partner,
        "multipartner_n_partners" = multipartners_n_partner,
        "prop_largest_multipartner" = sapply(multi_partner, largest_proportion),
        "prop_largest_nonsingleton" = sapply(non_single_names, largest_proportion),
        "prop_largest_all" = sapply(names(triple_list), largest_proportion),
        "prop_largest_multipartner_weighted" = weighted_proportion(multi_partner),
        "prop_largest_nonsingleton_weighted" = weighted_proportion(non_single_names),
        "prop_largest_all_weighted" = weighted_proportion(names(triple_list))
    )
}


options <- function() {

    #CLI args
    spec <- matrix(c(
        "help",  "h", "0", "logical",   "Print this help and exit",
        "input", "i", "1", "character", "Path to allele calls",
        "cores", "c", "1", "integer",   "Number of processor cores to use",
        "out",   "o", "1", "character", "Out path"
    ), byrow = TRUE, ncol = 5)

    opt <- getopt(spec)

    if (!is.null(opt$help)) {
        cat(getopt(spec, usage = TRUE))
        q(status = 1)
    }
    opt
}

process <- function(data, f, cores) {

    first <- 1:(ncol(data) - 2)

    triples <- mclapply(first, process_triples,
                        func = f, df = data, mc.cores = cores)

    names(triples) <- colnames(data)[first]

    stats <- data.frame(t(sapply(triples, extract_stats)))
    stats$left_gene <- names(triples)

    # reorder column names to put "left gene" in column 1
    reorder <- colnames(stats)[c(ncol(stats), 1:(ncol(stats) - 1))]
    stats <- stats[, reorder]

    stats

}

main <- function() {

    opt <- options()

    ###### Setup #####

    cgmlst <- read.table(opt$input, sep = ",", header = TRUE,
                         stringsAsFactors = FALSE, row.names = 1)

    ##### By Flank #####

    stats_by_flank <- process(cgmlst, get_triples_by_flank, opt$cores)

    ##### By Centre #####

    stats_by_centre <- process(cgmlst, get_triples_by_centre, opt$cores)

    write.csv(stats_by_flank,
              file = paste(opt$out, "triplet_flanks.csv", sep = "/"),
              row.names = FALSE, quote = FALSE)

    write.csv(stats_by_centre,
              file = paste(opt$out, "triplet_centres.csv", sep = "/"),
              row.names = FALSE, quote = FALSE)

}

main()
