library(randomForest)
library(magrittr)
library(getopt)

calculate_accuracy <- function(solution) {

    pick_result <- function(name) {

        sum(apply(solution[,c(name, 'actual')], 1,
              function(x) x[1] == x[2])) / nrow(solution)
    }

    lapply(colnames(solution)[colnames(solution) != "actual"], pick_result)
}

predict_centre <- function(chunk, name_list, train_size, trees) {

    set.seed(1) # for reproducibilty

    colnames(chunk) <- name_list

    selected <- sample(1:nrow(chunk), train_size, replace = FALSE)
    train <- chunk[selected, ]
    test <- chunk[-selected, ]

    solutions <- list()
    for (name in name_list) {
        if (name != 'centre')

            solutions[[name]] <-
                paste("as.factor(centre) ~", name) %>%
                as.formula %>%
                randomForest(data = train, ntree = trees) %>%
                predict(test)

    }

    if (length(name_list) > 3) {

        solutions[["lefts"]] <-
            paste(name_list[1:floor(ncol(chunk) / 2)], collapse = " + ") %>%
            paste("as.factor(centre) ~", .) %>%
            as.formula %>%
            randomForest(data = train, importance = TRUE, ntree = trees) %>%
            predict(test)


        solutions[["rights"]] <-
            name_list[(ceiling(ncol(chunk)/2) + 1):ncol(chunk)] %>%
            paste(collapse = " + ") %>%
            paste("as.factor(centre) ~", .) %>%
            as.formula %>%
            randomForest(data = train, importance = TRUE, ntree = trees) %>%
            predict(test)
    }

    solutions[["all"]] <-
        "as.factor(centre) ~" %>%
        paste(paste(name_list[name_list != 'centre'], collapse = " + ")) %>%
        as.formula %>%
        randomForest(data = train, importance = TRUE, ntree = trees) %>%
        predict(test)

    solutions[["actual"]] <- as.factor(test$centre)

    solution <- as.data.frame(solutions)

    solution


}

get_names <- function(width) {

    left <-  paste0("left_",  seq(floor(width / 2), 1))
    right <- paste0("right_", seq(1, floor(width / 2)))

    name_list <- c(left, "centre", right)

    name_list
}

generate_solutions <- function(calls, width, train, trees, cores) {

    generate_solution <- function(left) {

        chunk_indices <- left:(left + (width - 1))

        # ensure we aren't feeding in NAs
        chunk <- na.omit(calls[, chunk_indices])

        predict_centre(chunk = chunk,
                       name_list = get_names(width),
                       train_size = as.integer(nrow(calls) * train),
                       trees = trees)

    }

    1:(ncol(calls) - (width - 1)) %>%

    parallel::mclapply(generate_solution, mc.cores = cores )
}


load_calls <- function(input, min) {

    calls <-
        input %>%
        read.csv(stringsAsFactors = FALSE, row.names = 1)

    # drop columns with too few alleles
    pass_min <-
        calls %>%
        sapply(function(x) length(unique(x))) >= min

    calls %>% extract(pass_min)

}

options <- function() {
    # commandline options

    min_help <- "Minimum number of unique alleles for inclusion (default = 1)"
    width_help <- paste0("Width of region to be considered (must be odd);",
                         "Multiple values can be given as comma-seperated list")
    train_help <- "Proportion of data to be used as training (default = 0.1)"

    spec <- matrix(c(
        "help",  "h", "0", "logical",   "Print this help and exit",
        "input", "i", "1", "character", "Path to allele calls",
        "out",   "o", "1", "character", "Output path",
        "min",   "m", "2", "integer",   min_help,
        "width", "w", "1", "character", width_help,
        "train", "t", "2", "numeric",   train_help,
        "trees", "r", "1", "integer",   "Number of trees in the random forest",
        "cores", "c", "1", "integer",   "Number of CPU cores to use"
    ), byrow = TRUE, ncol = 5)

    opt <- getopt(spec)

    if (!is.null(opt$help)) {
        cat(getopt(spec, usage = TRUE))
        q(status = 1)
    }

    if (is.null(opt$minimum)) {
        opt$minimum <- 1
    }

    if (is.null(opt$train)) {
        opt$train <- 0.1
    }

    if (opt$width %% 2 == 0) {
        warn <- paste0("\n", "Argument --width must be odd!", "\n")
        cat(warn)
        q(status = 1)
    }

    opt$width <- as.integer(strsplit(opt$width, ",")[[1]])

    opt
}


main <- function() {

    opt <- options()

    calls <-
        opt$input %>%
        load_calls(opt$min)

    for (w in opt$width) {

        solutions <-
            calls %>%
            generate_solutions(w, opt$train, opt$trees, opt$cores)

        # percentage of the time the algo correctly predicted the allele
        assessments <-
            solutions %>%
            sapply(calculate_accuracy) %>%
            unlist %>%
            matrix(ncol = w + 2) %>%
            data.frame

        # fix data.frame headers
        sln1 <- solutions[[1]]
        colnames(assessments) <- colnames(sln1)[-ncol(sln1)]

        assessments %>% summary

        write.csv(assessments,
                  file = paste0(opt$out, "/width_", w, ".csv"),
                  row.names = FALSE, quote = FALSE)
    }
}
