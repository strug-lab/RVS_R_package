% Generated by roxygen2 (4.0.2.9000): do not edit by hand
\name{vp}
\alias{vp}
\title{Variance of genotypes}
\usage{
vp(P)
}
\arguments{
\item{P}{vector of 3 values for variant j, P(G=0), P(G=1) and P(G=2). (double)}
}
\value{
vp variance value (double)
}
\description{
This function computes the variance of a genotype given the 3 estimated genotype frequencies for variant j: P(G=0), P(G=1) and P(G=2).
}
\details{
vp is called from the \code{test_EGC} and \code{test_EGB} functions.
}
\seealso{
\code{\link{test_EGC}} and \code{\link{test_EGB}}
}

