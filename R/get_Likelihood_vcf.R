#' Gets Likelihoods
#'
#' This function fetches the genotype likelihoods from the VCF file (wo header) for the analysis and
#' must contain a 9th column (ie. 'FORMAT' field) consisting of (GT:AD:DP:GQ:PL).
#'
#' The first step is to read the vcf file as a matrix and outputs an R object with 5 lists:
#'
#' Input matrix:\cr
#' 2nd column (POS)    - location\cr
#' 4th column (REF)    - reference allele\cr
#' 5th column (ALT)    - alternative allele\cr
#' 7th column (FILTER) - vcr filter value\cr
#' 10th column and onwards are 'data' columns in format (GT:AD:DP:GQ:PL)
#'
#' PL = Phred-scaled likelihoods for genotypes as defined in the VCF specification format of 1000 Genome project.\cr
#' Phred-scaled likelihoods: log10(L(D|0)), log10(L(D|1)), log10(L(D|2)).
#'
#' Example: 0/0:2,0:4:6:0,6,42
#'
#' @param M The vcf file without headers in a matrix.
#' @return Outputs an R object with:\cr
#' 1. REF allele (vector)\cr
#' 2. ALT allele (vector)\cr
#' 3. POS of variant location (vector)\cr
#' 4. FILTER indicator for 'PASS' (vector)\cr
#' 5. likelihoods: L(D|0), L(D|1), L(D|2). List of 3 vectors.
#'
#' L0 - matrix of P(D|0), rows are individuals and columns are SNPS\cr
#' L1 - matrix of P(D|1), rows are individuals and columns are SNPS\cr
#' L2 - matrix of P(D|2), rows are individuals and columns are SNPS\cr
get_Likelihood_vcf <- function( M ){
	RA = NULL
	AA = NULL
	VL = NULL
	VI = NULL
	LD0 = NULL
	LD1 = NULL
	LD2 = NULL
	Lv = length(M[,1])
	Lind = length(M[1,])

	for (i in 1:Lv){
		RA = c(RA,as.vector(M[i,4]))
		AA = c(AA,as.vector(M[i,5]))
		VL = c(VL,M[i,2])
		VI = c(VI,as.vector(M[i,7]))
		l0 = NULL
		l1 = NULL
		l2 = NULL
		for (j in 10:Lind){
			if (M[i,j]=='./.') { LL = c(NA,NA,NA) }
			else { LL = get_l(M[i,j]) }
			l0 = c(l0,LL[1])
			l1 = c(l1,LL[2])
			l2 = c(l2,LL[3])
		}
		LD0 = cbind(LD0,l0)
		LD1 = cbind(LD1,l1)
		LD2 = cbind(LD2,l2)
	}
	return ( list( ref_a=RA, alt_a=AA, var_loc=VL, var_indic=VI, L_matrix=list( L0=LD0, L1=LD1, L2=LD2) ) )
}


#' Parse Phred-scaled likelihoods
#'
#' This function parses the data column (GT:AD:DP:GQ:PL) to get the PL values \code{-10log10(x)}, 
#' and converts them to a vector of likelihoods: L0,L1,L2. It's called within get_Likelihood_vcf().
#'
#' @param M matrix M[i,j] of matrix, ie. a single data unit of format: GT:AD:DP:GQ:PL.  Example: get_l('0/0:2,0:4:6:0,6,42')
#' 
#' @keywords get_Likelihood_vcf
#' @return a vector of three values: l0, l1, l2
#'
#' @seealso \code{\link{get_Likelihood_vcf}} for parent function.

get_l <- function( M ){
	a = unlist(strsplit(as.vector(M),':'))
	if (a[5]=='.') { return (c(NA,NA,NA)) }
	l = as.numeric(unlist(strsplit(a[5],',')))
	l0 = l[1]
	l1 = l[2]
	l2 = l[3]
	Pl0 = 10^(-l0/10)
	Pl1 = 10^(-l1/10)
	Pl2 = 10^(-l2/10)
	return ( c(Pl0,Pl1,Pl2) )
}


