% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PRS_PGx_Bayes.R
\name{PRS_PGx_Bayes}
\alias{PRS_PGx_Bayes}
\title{Construct PGx PRS using Bayesian regression}
\usage{
PRS_PGx_Bayes(
  PGx_GWAS,
  G_reference,
  n.itr = 1000,
  n.burnin = 500,
  n.gap = 10,
  paras,
  standardize = TRUE
)
}
\arguments{
\item{PGx_GWAS}{a numeric list containing PGx GWAS summary statistics (with SNP ID, position, \eqn{\beta}, \eqn{\alpha}, 2-df p-value, MAF and N), SD(Y), and mean(T)}

\item{G_reference}{a numeric matrix containing the individual-level genotype information from the reference panel (e.g., 1KG)}

\item{n.itr}{a numeric value indicating the total number of MCMC iteration}

\item{n.burnin}{a numeric value indicating the number of burn in}

\item{n.gap}{a numeric value indicating the MCMC gap}

\item{paras}{a numeric vector containg hyper-parameters (\eqn{v}, \eqn{\phi})}

\item{standardize}{a logical flag indicating should phenotype and genotype be standardized}
}
\value{
A numeric list, the first sublist contains estimated prognostic effect sizes, the second sublist contains estimated predictive effect sizes
}
\description{
Flexibly shrink prognostic and predictive effect sizes simutaneously with glocal-local shrinkage parameters
}
\details{
PRS-PGx-Bayes only needs PGx summary statistics and external reference genotype
}
\examples{
\donttest{
data(PRSPGx.example); attach(PRSPGx.example)
paras = c(3, 5)
coef_est <- PRS_PGx_Bayes(PGx_GWAS, G_reference, paras = paras, n.itr = 10, n.burnin = 5, n.gap = 1)
summary(coef_est$coef.G)
summary(coef_est$coef.TG)
}

}
\references{
Ge, T., Chen, CY., Ni, Y. et al. Polygenic prediction via Bayesian regression and continuous shrinkage priors. Nat. Commun. 10, 1776 (2019).

Zhai, S., Zhang, H., Mehrotra, D.V. & Shen, J. Paradigm Shift from Disease PRS to PGx PRS for Drug Response Prediction using PRS-PGx Methods (submitted).
}
\author{
Song Zhai
}
