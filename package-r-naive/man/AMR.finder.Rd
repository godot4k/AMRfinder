\name{AMRfinder}
\alias{AMRfinder}
\title{Perform AMR Segmentation Analysis}
\description{
  A function to perform segmentation analysis on DNA methylation data using the AMR method.
  This branch scores each candidate region with a phenotype-on-methylation model.
}
\usage{
  AMRfinder(intput_dat, y, cov.mod)
}
\arguments{
  \item{intput_dat}{A matrix containing input data with columns 'chr', 'pos', and CpG values.}
  \item{y}{A vector or data frame whose first column is the phenotype coded as 0/1. Additional
  columns are treated as sample-level adjustment variables, for example pre/post indicators.}
  \item{cov.mod}{An optional data frame of additional sample-level covariates to include in the model.}
  \item{controllist}{A list of control parameters}
}
\value{
  A data frame with segmentation results including chr, start, end, #CpGs, correlation estimates,
  methylation coefficients from phenotype ~ methylation + adjustments, p-values, methX,
  methY, metevalue-style e-value, FDR, BH-adjusted
  e-value, and e-BH significance flag.
}
\examples{
  # Example usage (if test data were available):
  # data <- readRDS("bulk.sub.txt.20.Rds")
  # result <- AMRfinder(data, y, cov.mod)
}
\author{Yifan Yang <yfyang.86@hotmail.com>}
\seealso{
  \code{\link[stats]{glm}} for phenotype model fitting.
}
\keyword{file}
