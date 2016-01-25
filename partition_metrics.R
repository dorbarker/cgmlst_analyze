
wallace <- function(a,b,c)
{
    # Compute the Wallace coefficients from a mismatch matrix from 2
    # classifications of a dataset.
    # Wallace, D. L. 1983. A method for comparing two hierarchical clusterings:
    # comment. J. Am. Stat. Assoc. 78:569-576.
    #
    # Args:
    #   a: a variable in formula to calculate the Wallace coefficient.
    #   b: b variable in formula to calculate the Wallace coefficient.
    #   c: c variable in formula to calculate the Wallace coefficient.
    #
    # Returns:
    #   Wallace coefficients for clustering 1 to 2 and 2 to 1.
    w1 = 0
    w2 = 0
    if ((a + b) > 0)
    {
        w1 = a / (a + b)
    }
    if ((a + c) > 0)
    {
        w2 = a / (a + c)
    }
    list(w1vs2 = w1, w2vs1 = w2)
}

simpsons <- function(ar)
{
    # Compute the Simpson's index of diversity for a list of classifications.
    # Simpson, E. H. 1949. Measurement of species diversity. Nature. 163:688.
    #
    # Args:
    #   ar: List of classifications.
    #
    # Returns:
    #   Simpson's index of diversity and 95% confidence interval values and the
    #   number of partitions.
    n <- length(ar)
    d <- table(ar)
    sid <- 1.0
    if (n * (n - 1) > 0)
    {
        sid <- 1.0 - (sum(d * (d - 1)) / (n * (n - 1)))
    }

    s2 <- 4 / n * (sum((d/n)^3) - sum((d/n)^2)^2)

    low <- sid - 2 * sqrt(s2)
    high <- sid + 2 * sqrt(s2)

    list(sid = sid,
         low = low,
         high = high,
         n_clusters = length(unique(ar))
        )
}

adj_wallace <- function(clusters_A, clusters_B)
{
    # Compute the Adjusted Wallace coefficients with 95% confidence intervals.
    # Severiano A., F. R. Pinto, M. Ramirez and J. Carriço. 2011. Adjusted
    #   Wallace as a Measure of Congruence between Typing Methods. J. Clin.
    #   Microbiol. doi:10.1128/JCM.00624-11
    # Pinto, F.R., J. Melo-Cristino, M. Ramirez. 2008. A Confidence Interval
    #   for the Wallace Coefficient of Concordance and Its Application to
    #   Microbial Typing Methods. PLoS ONE 3(11):e3696.
    #   doi:10.1371/journal.pone.0003696.
    #
    # Args:
    #   clusters_A: Array/list of one set of classifications.
    #   clusters_B: Array/list of another set of classifications.
    #
    # Returns:
    #   Adjusted Wallace coefficients with 95% confidence intervals for
    #   classifications 1 to 2 and 2 to 1.
    # print('adj_wallace')
    ct <- as.matrix(table(clusters_A, clusters_B))

    a <- sum(ct * (ct-1)/2)

    sum_col <- colSums(ct)

    sum_row <- rowSums(ct)

    a1 <- sum(sum_col * (sum_col -1 )/2)

    b <- a1 - a
    c <- sum(sum_row * (sum_row -1)/2) - a

    n <- sum(sum_row)
    d <- ((n * (n - 1)) / 2) - a1 - c


    sum_row_2 <- sum_row[sum_row > 1]
    sum_col_2 <- sum_col[sum_col > 1]
    csumFc2 <- colSums(sweep(ct[,names(sum_col_2)], 2, sum_col_2, FUN=function(x,y){(x/y)^2}))

    csumFc3 <- colSums(sweep(ct[,names(sum_col_2)], 2, sum_col_2, FUN=function(x,y){(x/y)^3}))

    rsumFc2 <- rowSums(sweep(ct[names(sum_row_2),], 1, sum_row_2, FUN=function(x,y){(x/y)^2}))

    rsumFc3 <- rowSums(sweep(ct[names(sum_row_2),], 1, sum_row_2, FUN=function(x,y){(x/y)^3}))

    rSumW1 <- sum((sum_row_2 * (sum_row_2 -1))^2 * (((4.0 * sum_row_2 * (sum_row_2 - 1.0) * (sum_row_2 - 2.0) * rsumFc3) + (2.0 * sum_row_2 * (sum_row_2 - 1.0) * rsumFc2) - (2.0 * sum_row_2 * (sum_row_2 - 1.0) * (2.0 * sum_row_2 - 3.0) * (rsumFc2^2.0))) / ((sum_row_2 * (sum_row_2 - 1.0))^2.0)))
    rSumW2 <- sum(sum_row * (sum_row - 1))
    cSumW1 <- sum((sum_col_2 * (sum_col_2 -1))^2 * (((4.0 * sum_col_2 * (sum_col_2 - 1.0) * (sum_col_2 - 2.0) * csumFc3) + (2.0 * sum_col_2 * (sum_col_2 - 1.0) * csumFc2) - (2.0 * sum_col_2 * (sum_col_2 - 1.0) * (2.0 * sum_col_2 - 3.0) * (csumFc2^2.0))) / ((sum_col_2 * (sum_col_2 - 1.0))^2.0)))
    cSumW2 <- sum(sum_col * (sum_col - 1))

    # get variance of Wallace 1vs2 (varW1) and 2vs1 (varW2)

    varW1 <- 0.0
    varW2 <- 0.0
    if (rSumW2 > 0)
    {
        varW1 <- (rSumW1 / ((rSumW2) ^ 2.0))
    }
    if (cSumW2 > 0)
    {
        varW2 <- (cSumW1 / ((cSumW2) ^ 2.0))
    }
    # get mismatch matrix
    wallace <- wallace(a, b, c)
    w1 <- wallace$w2vs1
    w2 <- wallace$w1vs2
    sid1 <- simpsons(clusters_A)
    sid2 <- simpsons(clusters_B)

    wi1 <- 1 - sid1$sid
    wi2 <- 1 - sid2$sid
    aw1 <- (w1 - wi2) / (1 - wi2)
    aw2 <- (w2 - wi1) / (1 - wi1)
    aw1CI <- 2.0 * (1.0 / (1.0 - wi2)) * sqrt(varW1)
    aw1Low <- aw1 - aw1CI
    if (aw1Low < 0.0)
    {
        aw1Low <- 0.0
    }
    aw1High <- aw1 + aw1CI
    if (aw1High > 1.0)
    {
        aw1High <- 1.0
    }
    aw2CI <- 2.0 * (1.0 / (1.0 - wi1)) * sqrt(varW2)
    aw2Low <- aw2 - aw2CI
    aw2High <- aw2 + aw2CI
    if (aw2Low < 0.0)
    {
        aw2Low <- 0.0
    }
    if (aw2High > 1.0)
    {
        aw2High <- 1.0
    }

    list(Wallace_A_vs_B = w1,
         Wallace_B_vs_A = w2,
         Simpsons_A = sid1,
         Simpsons_B = sid2,
         Adjusted_Wallace_A_vs_B = aw1,
         Adjusted_Wallace_A_vs_B_low = aw1Low,
         Adjusted_Wallace_A_vs_B_high = aw1High,
         Adjusted_Wallace_B_vs_A = aw2,
         Adjusted_Wallace_B_vs_A_low = aw2Low,
         Adjusted_Wallace_B_vs_A_high = aw2High
        )

}

taboada_cohesion <- function(clustersA, clustersB) {

    mat <- table(clustersA, clustersB)

    strains_in_clusts <- apply(mat, 1, sum)

    size_of_largest_subclust <- apply(mat, 1, max)

    corresponding_clust <- unlist(apply(mat, 1, function(x) which(x == max(x))[1]))

    proportion <- apply(cbind(strains_in_clusts, size_of_largest_subclust), 1, function(x) x[2] / x[1])

    concordant <- apply(cbind(strains_in_clusts, proportion), 1, function(x) x[1] * x[2])

    cluster_hits <- apply(mat, 1, function(x) sum(x > 0))

    weighted_avg <- sum(concordant) / sum(strains_in_clusts)

    list("n_strains" = strains_in_clusts,
         "largest_size" = size_of_largest_subclust,
         "corresponding" = corresponding_clust,
         "proportions" = proportion,
         "n_concordant" = concordant,
         "cluster_hits" = cluster_hits,
         "weighted_average" = weighted_avg)

}
