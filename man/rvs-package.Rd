\name{rvs-package}
\alias{rvs-package}
\alias{rvs}
\docType{package}
\title{
	Association Analysis using Next-Gen Sequence Data.
}
\description{
	The robust variance score (RVS) is a likelihood-based method/test to compute
	p-values to analyze and interpret genetic association data.  It
	substitutes genotype calls by their expected values given observed sequence
	data.

	RVS eliminates read depth bias in the estimation of minor allele
	frequency. It controls Type I error and has comparable power to the 'gold
	standard' analysis with the true underlying genotypes for both common and
	rare variants.

	The current package contains all the core functions but does not process or
	filter the vcf file yet.

}
\details{
	\tabular{ll}{
		Package: \tab rvs\cr
		Type: \tab Package\cr
		Version: \tab 0.0\cr
		Date: \tab 2014-10-20\cr
		License: \tab GPL (>= 3)\cr
	}
	The main Regular Score/Trend test functions are \code{test_EG} and \code{test_EGPv}
	which by evaluated asymptotic and permutation distributions respectively. The main RVS
	test functions are \code{test_EGC}, \code{test_EGC}, \code{test_EGC}

}
\author{
	Dr. Lisa J Strug <lisa.strug@utoronto.ca>\cr
	Andriy Derkach <derkach@utstat.toronto.edu>\cr
	Maintainer: Ted Chiang <tchiang@sickkids.ca>\cr
}
\references{
Derkach et. al, Bioinformatics: Vol.30 no.15 2014, pp.2179 to 2188.
}
\seealso{
 \code{\link{test_EG}},
 \code{\link{test_EGPv}},
 \code{\link{test_EGC}},
 \code{\link{test_EGCV}}
 \code{\link{test_EGB}},
 \code{\link{test_EGBnew}},
}
\keyword{models}

