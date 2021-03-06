% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_for_RGWAS.R
\name{score.calc}
\alias{score.calc}
\title{Calculate -log10(p) for single-SNP GWAS}
\usage{
score.calc(
  M.now,
  ZETA.now,
  y,
  X.now,
  Hinv,
  P3D = TRUE,
  optimizer = "nlminb",
  eigen.G = NULL,
  n.core = NA,
  min.MAF = 0.02,
  count = TRUE
)
}
\arguments{
\item{M.now}{A n.sample x n.mark genotype matrix where n.sample is sample size and n.mark is the number of markers.}

\item{ZETA.now}{A list of variance (relationship) matrix (K; \eqn{m \times m}) and its design matrix (Z; \eqn{n \times m}) of random effects. You can use only one kernel matrix.
For example, ZETA = list(A = list(Z = Z, K = K))
Please set names of list "Z" and "K"!}

\item{y}{A \eqn{n \times 1} vector. A vector of phenotypic values should be used. NA is allowed.}

\item{X.now}{A \eqn{n \times p} matrix. You should assign mean vector (rep(1, n)) and covariates. NA is not allowed.}

\item{Hinv}{The inverse of \eqn{H = ZKZ' + \lambda I} where \eqn{\lambda = \sigma^2_e / \sigma^2_u}.}

\item{P3D}{When P3D = TRUE, variance components are estimated by REML only once, without any markers in the model.
When P3D = FALSE, variance components are estimated by REML for each marker separately.}

\item{optimizer}{The function used in the optimization process. We offer "optim", "optimx", and "nlminb" functions.}

\item{eigen.G}{A list with
\describe{
\item{$values}{Eigen values}
\item{$vectors}{Eigen vectors}
}
The result of the eigen decompsition of \eqn{G = ZKZ'}. You can use "spectralG.cpp" function in RAINBOWR.
If this argument is NULL, the eigen decomposition will be performed in this function.
We recommend you assign the result of the eigen decomposition beforehand for time saving.}

\item{n.core}{Setting n.core > 1 will enable parallel execution on a machine with multiple cores.}

\item{min.MAF}{Specifies the minimum minor allele frequency (MAF).
If a marker has a MAF less than min.MAF, it is assigned a zero score.}

\item{count}{When count is TRUE, you can know how far RGWAS has ended with percent display.}
}
\value{
-log10(p) for each marker
}
\description{
Calculate -log10(p) of each SNP by the Wald test.
}
\references{
Kennedy, B.W., Quinton, M. and van Arendonk, J.A. (1992)
Estimation of effects of single genes on quantitative traits. J Anim Sci. 70(7): 2000-2012.

Kang, H.M. et al. (2008) Efficient Control of Population Structure
 in Model Organism Association Mapping. Genetics. 178(3): 1709-1723.

Kang, H.M. et al. (2010) Variance component model to account for sample
  structure in genome-wide association studies. Nat Genet. 42(4): 348-354.

Zhang, Z. et al. (2010) Mixed linear model approach adapted for genome-wide
 association studies. Nat Genet. 42(4): 355-360.
}
