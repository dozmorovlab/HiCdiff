% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hic_loess.R
\name{hic_loess}
\alias{hic_loess}
\title{Perform joint loess normalization on two Hi-C datasets}
\usage{
hic_loess(hic.table, degree = 1, span = NA, loess.criterion = "gcv",
  Plot = FALSE, parallel = FALSE, BP_param = bpparam(),
  check.differences = FALSE, diff.thresh = "auto", iterations = 10000)
}
\arguments{
\item{hic.table}{hic.table or a list of hic.tables generated from
 the create.hic.table function.
 list of hic.tables generated from the create.hic.table function.
 If you want to perform
normalization over multiple chromosomes from each cell line at
once utilizing parallel computing enter a list of hic.tables and
set parallel = TRUE.}

\item{degree}{Degree of polynomial to be used for loess. Options
are 0, 1, 2. The default setting is degree = 1.}

\item{span}{User set span for loess. If set to NA, the span will
be selected automatically using the setting of loess.criterion.
Defaults to NA so that
automatic span selection is performed.}

\item{loess.criterion}{Automatic span selection criterion. Can use either
'gcv' for generalized cross-validation or 'aicc' for Akaike Information
 Criterion.
 Span selection uses a slightly modified version of the \code{loess.as()}
 function from the \code{fANCOVA} package. Defaults to 'gcv'.}

\item{Plot}{Logical, should the MD plot showing before/after loess
normalization be output? Defaults to FALSE.}

\item{parallel}{Logical, set to TRUE to utilize the \code{parallel} package's
parallelized computing. Only works on unix operating systems. Only useful if
entering a list of hic.tables. Defauts to FALSE.}

\item{BP_param}{Parameters for BiocParallel. Defaults to bpparam(), see help
for BiocParallel for more information
\url{http://bioconductor.org/packages/release/bioc/vignettes/BiocParallel/
inst/doc/Introduction_To_BiocParallel.pdf}}

\item{check.differences}{Logical, should difference detection be performed? If TRUE,
the same procedure as hic_diff will be performed. If FALSE,
only normalization will be performed on the entered data. Defaults to FALSE.}

\item{diff.thresh}{Fold change threshold desired to call a detected difference
clinically significant. Set to "auto" by default to indicate that the
difference threshold will be automatically calculated as 2 standard deviations
of all the adjusted M values. For no p-value adjustment
set diff.thresh = NA. To set your own threshold enter a numeric value i.e.
diff.thresh = 1. If set to "auto" or a numeric value, a check will
be made as follows: if permutation p-value < 0.05 AND M < diff.thresh (the log2
fold change for the difference between IF1 and IF2) then
the p-value will be set to 0.5. Defaults to 'auto'.}

\item{iterations}{Number of iterations for the permuation test. Will only be used
if check.differences set to TRUE. Defaults to 10000.}
}
\value{
An updated hic.table is returned with the additional columns of adj.IF1,
   adj.IF2 for the respective normalized IFs, an adj.M column for the
   adjusted M, and mc for the loess correction factor. If \code{check.differences}
   is set to TRUE a column containing the p-values for the
   significance of the difference between the two datasets will also be returned.
}
\description{
Perform joint loess normalization on two Hi-C datasets
}
\details{
The function takes in a hic.table or a list of hic.table objects created
    with the \code{create.hic.loess} function. If you wish to perform joint
    normalization on Hi-C data for multiple chromosomes use a list of hic.tables.
    The process can be parallelized using the \code{parallel}
    setting. The data is fist transformed into what is termed an MD plot (similar
    to the MA plot/Bland-Altman plot). M is the log difference log2(x/y) between
    the two datasets. D is the unit distance in the contact matrix. The MD plot can
    be visualized with the \code{Plot} option. Loess regression is then
    performed on the MD plot to model any biases between the two Hi-C datasets. An
    adjusted IF is then calculated for each dataset along with an adjusted M.
    See methods section of Stansfield & Dozmorov 2017 for more details. Note:
    if you receive the warning "In simpleLoess(y, x, w, span, degree = degree,
    parametric = parametric,  ... :pseudoinverse used..." it should not effect
    your results, however it can be avoided by manually setting the span to
    a larger value using the span option.
}
\examples{
# Create hic.table object using included Hi-C data in sparse upper
# triangular matrix format
data("HMEC.chr22")
data("NHEK.chr22")
hic.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr= 'chr22')
# Plug hic.table into hic_loess()
result <- hic_loess(hic.table, Plot = TRUE)
# View result
result

}
