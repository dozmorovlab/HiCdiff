#' Detect differences between two jointly normalized Hi-C datasets.
#'
#' @param hic.table A hic.table or list of hic.tables output from the
#'     \code{hic_loess} function. hic.table must be jointly normalized
#'     before being entered.
#' @param diff.thresh Fold change threshold desired to call a detected
#'     difference clinically significant. Set to 'auto' by default to indicate that the
#'     difference threshold will be automatically calculated as 2 standard
#'     deviations of all the adjusted M values. For no p-value adjustment
#'     set diff.thresh = NA. To set your own threshold enter a numeric value
#'     i.e. diff.thresh = 1. If set to 'auto' or a numeric value, a check will
#'     be made as follows: if permutation p-value < 0.05 AND M < diff.thresh (the log2
#'     fold change for the difference between IF1 and IF2) then
#'     the p-value will be set to 0.5.
#' @param iterations Number of iterations for the permuation test.
#' @param Plot Logical, should the MD plot showing before/after loess normalization
#'     be output?
#' @param parallel Logical, set to TRUE to utilize the \code{parallel} package's
#'     parallelized computing. Only works on unix operating systems. Only useful if
#'     entering a list of hic.tables.
#' @param numCores Number of cores to be used if parallel set to TRUE.
#'
#'
#' @details  The function takes in a hic.table or a list of hic.table objects created
#'     with the \code{hic_loess} function. If you wish to perform difference
#'     detection on Hi-C data for multiple chromosomes use a list of hic.tables. The process
#'     can be parallelized on unix systems using the \code{parallel}
#'     setting. The adjusted IF and adjusted M calculated from \code{hic_loess} are used for
#'     difference detection. A permutation test is performed to test
#'     the significance of the difference between each IF of the two datasets. Permutations
#'     are broken in blocks for each unit distance. See methods section
#'     of Stansfield & Dozmorov 2017 for more details.
#'
#' @return A hic.table with additional columns containing a p-value for the significance
#'     of the difference and the raw fold change between the IFs of the two datasets.
#'
#' @examples
#' # Create hic.table object using included Hi-C data in sparse upper triangular
#' # matrix format
#' data('HMEC.chr22')
#' data('NHEK.chr22')
#' hic.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr = 'chr22')
#' # Plug hic.table into hic_loess()
#' result <- hic_loess(hic.table, Plot = TRUE)
#' # perform difference detection
#' diff.result <- hic_diff(result, diff.thresh = 'auto', Plot = TRUE)
#'
hic_diff <- function(hic.table, diff.thresh = "auto", iterations = 10000,
                     Plot = FALSE, parallel = FALSE, numCores = 3) {
  # check for correct input
  if (!is.na(diff.thresh) & is.numeric(diff.thresh) & diff.thresh <=
      0) {
    stop("Enter a numeric value > 0 for diff.thresh or set it to NA or \"auto\"")
  }
  if (!is.na(diff.thresh) & is.character(diff.thresh) & diff.thresh !=
      "auto") {
    stop("Enter a numeric value > 0 for diff.thresh or set it to NA or \"auto\"")
  }
  # check if single hic.table or list
  if (is.data.table(hic.table)) {
    hic.table <- list(hic.table)
  }
  # calculate diff.thresh if set to auto
  if (!is.na(diff.thresh)) {
    if (diff.thresh == "auto") {
      diff.thresh <- sapply(hic.table, .calc.diff.thresh)
    }
  }
  # run difference detection for parallel / non-parallel
  if (parallel) {
    if (length(diff.thresh) == 1) {
      hic.table <- mclapply(hic.table, .calc.pval, Plot = Plot, diff.thresh = diff.thresh,
                            iterations = iterations)
    } else {
      hic.table <- mcmapply(.calc.pval, hic.table, diff.thresh,
                            MoreArgs = list(Plot = Plot,
                                            iterations = iterations), SIMPLIFY = FALSE)
    }
  } else {
    if (length(diff.thresh) == 1) {
      hic.table <- lapply(hic.table, .calc.pval, Plot = Plot,
                          diff.thresh = diff.thresh,
                          iterations = iterations)
    } else {
      hic.table <- mapply(.calc.pval, hic.table, diff.thresh,
                          MoreArgs = list(Plot = Plot,
                                          iterations = iterations), SIMPLIFY = FALSE)
    }
  }
  # clean up if single hic.table
  if (length(hic.table) == 1) {
    hic.table <- hic.table[[1]]
  }
  return(hic.table)
}
