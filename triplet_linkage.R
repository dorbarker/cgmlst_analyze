library(parallel)

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
    for(i in 1:nrow(df)) {

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
# Setup
cgmlst <- read.table("~/Dropbox/2014_05_to_2014_08_Dillon_wgMLST_approaches/cgmlst_paper/allele_calls/cgmlst_calls.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
rownames(cgmlst) <- cgmlst[,1]
cgmlst <- cgmlst[,-1]
first_col <- 1:(ncol(cgmlst)-2)

## By Flank

flank_triples <- mclapply(first_col, process_triples, func = get_triples_by_flank, df = cgmlst, mc.cores = 23)
names(flank_triples) <- colnames(cgmlst)[first_col]

stats_by_flank <- data.frame(t(sapply(flank_triples, extract_stats)))
stats_by_flank$left_gene <- names(flank_triples)
stats_by_flank <- stats_by_flank[,colnames(stats_by_flank)[c(ncol(stats_by_flank), 1:(ncol(stats_by_flank)-1))]]

# By Centre

centre_triples <- mclapply(first_col, process_triples, func = get_triples_by_centre, df = cgmlst, mc.cores = 23)
names(centre_triples) <- colnames(cgmlst)[first_col]

stats_by_centre <- data.frame(t(sapply(centre_triples, extract_stats)))
stats_by_centre$left_gene <- names(centre_triples)
stats_by_centre <- stats_by_centre[,colnames(stats_by_centre)[c(ncol(stats_by_centre), 1:(ncol(stats_by_centre)-1))]]

