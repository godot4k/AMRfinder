\name{AMRfinder}
\alias{AMRfinder}
\title{Perform AMR Segmentation Analysis}
\description{
  A function to perform segmentation analysis on DNA methylation data using the AMR method.
  This branch scores each candidate region with a continuous phenotype-on-methylation
  linear or linear mixed model.
}
\usage{
  AMRfinder(intput_dat, y, cov.mod)
}
\arguments{
  \item{intput_dat}{A matrix containing input data with columns 'chr', 'pos', and CpG values.}
  \item{y}{A vector or data frame whose first column is a numeric continuous phenotype. Additional
  columns are treated as sample-level adjustment variables, for example pre/post indicators or
  subject identifiers.}
  \item{cov.mod}{An optional data frame of additional sample-level covariates to include in the model.}
  \item{controllist}{A list of control parameters}
}
\value{
  A data frame with segmentation results including chr, start, end, #CpGs, correlation estimates,
  methylation coefficients from phenotype ~ methylation + adjustments, p-values, methX,
  methY, and FDR.
}
\examples{
  # Example usage (if test data were available):
  # data <- readRDS("bulk.sub.txt.20.Rds")
  # result <- AMRfinder(data, y, cov.mod)
}
\author{Yifan Yang <yfyang.86@hotmail.com>}
\seealso{
  \code{\link[stats]{lm}} and \code{\link[lmerTest]{lmer}} for phenotype model fitting.
}
\keyword{file}
