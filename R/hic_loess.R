#' Perform joint loess normalization on two Hi-C datasets
#'
#' @param hic.table hic.table or a list of hic.tables generated from
#'     the create.hic.table function.
#'     list of hic.tables generated from the create.hic.table function.
#'     If you want to perform
#'    normalization over multiple chromosomes from each cell line at
#'    once utilizing parallel computing enter a list of hic.tables and
#'    set parallel = TRUE.
#' @param degree Degree of polynomial to be used for loess. Options
#'     are 0, 1, 2. The default setting is degree = 1.
#' @param span User set span for loess. If set to NA, the span will
#'     be selected automatically using the setting of loess.criterion.
#'     Defaults to NA so that
#'     automatic span selection is performed.
#' @param loess.criterion Automatic span selection criterion. Can use either
#'     'gcv' for generalized cross-validation or 'aicc' for Akaike Information
#'      Criterion.
#'      Span selection uses a slightly modified version of the \code{loess.as()}
#'      function from the \code{fANCOVA} package. Defaults to 'gcv'.
#' @param Plot Logical, should the MD plot showing before/after loess
#'     normalization be output? Defaults to FALSE.
#' @param parallel Logical, set to TRUE to utilize the \code{parallel} package's
#'     parallelized computing. Only works on unix operating systems. Only useful if
#'     entering a list of hic.tables. Defauts to FALSE.
#' @param numCores Number of cores to be used if parallel set to TRUE. Defaults to 3.
#' @param check.differences Logical, should difference detection be performed? If TRUE,
#'     the same procedure as hic_diff will be performed. If FALSE,
#'     only normalization will be performed on the entered data. Defaults to FALSE.
#' @param diff.thresh Fold change threshold desired to call a detected difference
#'     clinically significant. Set to "auto" by default to indicate that the
#'     difference threshold will be automatically calculated as 2 standard deviations
#'     of all the adjusted M values. For no p-value adjustment
#'     set diff.thresh = NA. To set your own threshold enter a numeric value i.e.
#'     diff.thresh = 1. If set to "auto" or a numeric value, a check will
#'     be made as follows: if permutation p-value < 0.05 AND M < diff.thresh (the log2
#'     fold change for the difference between IF1 and IF2) then
#'     the p-value will be set to 0.5. Defaults to 'auto'.
#' @param iterations Number of iterations for the permuation test. Will only be used
#'     if check.differences set to TRUE. Defaults to 10000.
#'
#' @details The function takes in a hic.table or a list of hic.table objects created
#'     with the \code{create.hic.loess} function. If you wish to perform joint
#'     normalization on Hi-C data for multiple chromosomes use a list of hic.tables.
#'     The process can be parallelized on unix systems using the \code{parallel}
#'     setting. The data is fist transformed into what is termed an MD plot (similar
#'     to the MA plot/Bland-Altman plot). M is the log difference log2(x/y) between
#'     the two datasets. D is the unit distance in the contact matrix. The MD plot can
#'     be visualized with the \code{Plot} option. Loess regression is then
#'     performed on the MD plot to model any biases between the two Hi-C datasets. An
#'     adjusted IF is then calculated for each dataset along with an adjusted M.
#'     See methods section of Stansfield & Dozmorov 2017 for more details. Note:
#'     if you receive the warning "In simpleLoess(y, x, w, span, degree = degree,
#'     parametric = parametric,  ... :pseudoinverse used..." it should not effect
#'     your results, however it can be avoided by manually setting the span to
#'     a larger value using the span option.
#'
#'
#' @return An updated hic.table is returned with the additional columns of adj.IF1,
#'    adj.IF2 for the respective normalized IFs, an adj.M column for the
#'    adjusted M, and mc for the loess correction factor. If \code{check.differences}
#'    is set to TRUE a column containing the p-values for the
#'    significance of the difference between the two datasets will also be returned.
#'
#' @examples
#' # Create hic.table object using included Hi-C data in sparse upper
#' # triangular matrix format
#' data("HMEC.chr22")
#' data("NHEK.chr22")
#' hic.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr= 'chr22')
#' # Plug hic.table into hic_loess()
#' result <- hic_loess(hic.table, Plot = TRUE)
#' # View result
#' result
#'
hic_loess = function(hic.table, degree = 1, span = NA, loess.criterion = 'gcv',
                     Plot = FALSE, parallel = FALSE, numCores = 3,
                     check.differences = FALSE, diff.thresh = 'auto',  iterations = 10000) {
  # check for correct inputs
  if (!is.na(diff.thresh) & is.numeric(diff.thresh) & diff.thresh <= 0) {
    stop('Enter a numeric value > 0 for diff.thresh or set it to NA or "auto"')
  }
  if (!is.na(diff.thresh) & is.character(diff.thresh) & diff.thresh != "auto") {
    stop('Enter a numeric value > 0 for diff.thresh or set it to NA or "auto"')
  }
  if (!is.na(span) & span < 0.001) {
    stop('Enter a larger value for span')
  }
  if (!is.na(span) & span > 1) {
    stop('Enter a value <= 1 for span')
  }
  # check if list or single hic.table
  if (is.data.table(hic.table)) hic.table = list(hic.table)
  # apply loess.matrix to the list of matrices
  if (parallel) hic.table = mclapply(hic.table, .loess.matrix, mc.cores = numCores,
                                     Plot = Plot, degree = degree, span = span,
                                     loess.criterion = loess.criterion)
  else hic.table = lapply(hic.table, .loess.matrix, Plot = Plot, degree = degree,
                          span = span, loess.criterion = loess.criterion)

  # calculate p-value matrices if option is TRUE
  if(check.differences) {
    if (!is.na(diff.thresh)) {
      if (diff.thresh == 'auto') {
        diff.thresh = sapply(hic.table, .calc.diff.thresh)
      }
    }
    if (parallel) {
      if (length(diff.thresh) == 1) {
        hic.table = mclapply(hic.table, .calc.pval, Plot = Plot,
                             diff.thresh = diff.thresh,
                             iterations = iterations)
      } else {
        hic.table = mcmapply(.calc.pval, hic.table, diff.thresh,
                             MoreArgs = list(Plot = Plot,
                                             iterations = iterations), SIMPLIFY = FALSE)
      }
    } else {
      if (length(diff.thresh) == 1) {
        hic.table = lapply(hic.table, .calc.pval, Plot = Plot,
                           diff.thresh = diff.thresh,
                           iterations = iterations)
      } else {
        hic.table = mapply(.calc.pval, hic.table, diff.thresh,
                           MoreArgs = list(Plot = Plot,
                                           iterations = iterations), SIMPLIFY = FALSE)
      }
    }
  }
  # if there is only one hic.table entered pull it out the list before returning it
  if (length(hic.table) == 1) hic.table = hic.table[[1]]
  return(hic.table)
}
