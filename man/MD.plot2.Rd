% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{MD.plot2}
\alias{MD.plot2}
\title{Visualize the MD plot.}
\usage{
MD.plot2(M, D, p.val = NA, diff.thresh = NA)
}
\arguments{
\item{M}{The M component of the MD plot. Available from the
hic.table object.}

\item{D}{The D component of the MD plot. The unit distance of the
interaction. Available from the hic.table object.}

\item{p.val}{An optional p-value vector to provide color to the plot
based on the significance of the differences between the IFs.}

\item{diff.thresh}{A difference threshold used for calculating p-values.
If set to a value will add dotted horizontal lines to the plot to
display the threshold cutoffs. See `hic_loess` or `hic_diff`
functions help for more information on this parameter.}
}
\value{
An MD plot.
}
\description{
Visualize the MD plot.
}
\examples{
data('HMEC.chr22')
data('NHEK.chr22')
hic.table <- create.hic.table(HMEC.chr22, NHEK.chr22, chr='chr22')
# Plug hic.table into hic_loess()
result <- hic_loess(hic.table)
# perform difference detection
diff.result <- hic_diff(result, diff.thresh = 1)
MD.plot2(diff.result$M, diff.result$D, diff.result$p.value)

}
