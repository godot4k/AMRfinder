\name{AMRfinder}
\alias{AMRfinder}
\title{Perform AMR Segmentation Analysis}
\description{
  A function to perform segmentation analysis on DNA methylation data using the AMR method.
}
\usage{
  AMRfinder(intput_dat, y, cov.mod = NULL, controlist = list())
}
\arguments{
  \item{intput_dat}{A matrix containing input data with columns 'chr', 'pos', and CpG values.}
  \item{y}{A vector or matrix representing the binary response variable. Control samples must be coded as 0 and test samples must be coded as 1.}
  \item{cov.mod}{Unused in the no.cov control/test version.}
  \item{controlist}{A list of control parameters.}
}
\value{
  A data frame with segmentation results including chr, start, end, N.CpGs,
  ks_stat, mean_diff, p_value, methX, methY, e_value, FDR, e_adjust,
  and e_bh_significant. The e_value column uses the metevalue-style
  group-specific versus pooled normal density ratio. The e_adjust column is
  calculated as 1 / p.adjust(1 / e_value, method = "BH").
}
\examples{
  # Example usage (if test data were available):
  # data <- readRDS("bulk.sub.txt.20.Rds")
  # y <- data.frame(y = c(rep(0, 8), rep(1, 8)))
  # result <- AMRfinder(data, y, NULL)
}
\author{Yifan Yang <yfyang.86@hotmail.com>}
\seealso{
  \code{\link[stats]{ks.test}} for the two-sample Kolmogorov-Smirnov test.
}
\keyword{file}
