#' Produces base calls
#'
#' This function produces the base call from given a genotype and the error scores. It is called from the \code{generate_seqdata} functions.
#'
#' @param genotype vector of two alleles for a given genotype, ie. c('C', 'T')
#' @param error vector of base call errors: Phen scores (double).  Length of this vector is the read depth.
#' @param gen_vect vector of prior genotypes.
#' @return vect_row_reads vector of row sequence reads ('A','G','T','C') for single variant. The length of this vector is the read depth.
#' @seealso \code{\link{get_row_reads}} 
seq_call <- function( genotype, error, gen_vect ){
	depth = length(error)
	vect_row_reads = get_row_reads( error, genotype, depth )
	return (vect_row_reads)
}

#' @describeIn seq_call see \code{seq_call} is the wrapper around get_row_reads.
#' @param depth the read depth (integer)
#' @seealso \code{\link{seq_call}} 
get_row_reads <- function( error, genotype, depth ){
	vect_row_reads = NULL
	etal = c('A','G','T','C')
	for (i in 1:depth){
		s = rbinom(1,1,0.5)+1
		g = genotype[s]
		ng = etal[etal!=g]
		e3 = error[i]/3
		e23 = e3*2
		e = e3*3
		r = runif(1)
		if (r>e) {vect_row_reads = c(vect_row_reads,g)}
		if (r<e3){vect_row_reads = c(vect_row_reads,ng[1])}
		if (r>e3 & r<e23){vect_row_reads = c(vect_row_reads,ng[2])}
		if (r>e23 & r<e){vect_row_reads = c(vect_row_reads,ng[3])}
	}
	return (vect_row_reads)
}




#' Get likelihood P(Di|G) matrix
#'
#' This function calculates a likelihood matrix given a vector of error rates and a vector of sequence reads.
#'
#' Called inside the \code{generate_seqdata} functions.  The output of this function \code{get_M} are used as inputs to \code{get_EG}.
#' Note to Andriy:  Why is LG hardcoded to LG=c('TT','CT','CC')  ??  Guess it's to denote the AA, Aa, aa.
#'
#' @param Error vector of base call error rates; Uses the same \code{error} input in the \code{seq_call} function but concatenate all \code{error} values over all samples.
#' @param vect_row_reads vector of row sequence reads ('A','G','T','C') outputted from function \code{seq_call}, but concatenated over all samples in savme fashion as \code{Error}.
#' @return M matrix of P(read|AA), P(read|Aa) and P(read|aa), where row dimensions equal \code{length(Error)}, by 3 columns.
#' @seealso \code{\link{seq_call}} 
get_M <- function( Error, vect_row_reads ){
	L = length(Error)
	M = NULL
	LG  =c('TT','CT','CC')

	for (i in 1:L){
		m = NULL
		for  (j in 1:3){
			G = LG[j]
			A1 = substring(G,1,1)
			A2 = substring(G,2,2)
			m = c(m,get_Lsingle(vect_row_reads[i],A1,A2,Error[i]))
		}		
		M = rbind(M,m)
	}
	return (M)
}

#' @describeIn get_M with \code{rdv}
#' @param rdv vector of read depths concatenated for all samples (integer).   Looks like getting likelihood over all reads cumulatively? Andriy, please explain.
get_Mr <- function( Error, vect_row_reads, rdv ){
	L = length(rdv)
	M = NULL
	LG  =c('TT','CT','CC')

	S = 0
	for (i in 1:L){
		m = NULL
		for  (j in 1:3){
			LL = 1
			for (kk in 1:rdv[i]){
				G = LG[j]
				A1 = substring(G,1,1)
				A2 = substring(G,2,2)
				LL = LL*get_Lsingle(vect_row_reads[S+kk],A1,A2,Error[S+kk])
			}
			m =c(m,LL)
		}
		S =S+rdv[i]
		M = rbind(M,m)
	}
	return (M)
}

#' Get the genotype likelihood for a single sequence read
#'
#' This function computes the genotype likelihood for a single base read given the sequencing error rate.
#' 
#' get_Lsingle is called from the \code{get_M} and \code{get_Mr} functions.
#'
#' @param vect_row_read single read, ie. ('A')
#' @param A1 reference allele, ie. (A1-'A',A2-'C')
#' @param A2 reference allele, ie. (A1-'A',A2-'C')
#' @param error sequencing error (double)
#' @return p - likelihood of A1 and A2 for single read
#' @seealso \code{\link{get_M}} or \code{\link{get_Mr}} 
get_Lsingle <- function( vect_row_read, A1, A2, error ){
	if (vect_row_read == A1){ p1 = 1-error }
	else{ p1 = error/3 }

	if (vect_row_read == A2){ p2 = 1-error }
	else{ p2 = error/3 }

	p = 1/2*p1 + 1/2*p2
	return (p)
}




#' Conditional Expected value of genotype
#'
#' This function calculates the conditional expected value (or conditional variance value) of the genotype given the genotype likelihoods and frequencies.
#'
#' @param M genotype likelihoods AA, Aa, aa, matrix sum(rdv) by 3 (double); uses output from get_M function.
#' @param p genotype frequencies AA, Aa, aa (double); output from EM_calc function.
#' @param rdv read depth (vector of integers) for all samples
#'
#' @return EG a vector containing conditional expectation values. Size should be equal to length(rdv).
get_EG <- function( M, p, rdv ){
	L = length(rdv)
	S = 0
	EG = NULL
	g = c(0,1,2)
	for (i in 1:L){
		m = NULL
		for  (j in 1:3){
			L = 1
			for (kk in 1:rdv[i]){
				L = L*M[S + kk,j]
			}
			m = c(m,L*p[j])
		}
		S = S + rdv[i]
		pm = sum(m/sum(m)*g)
		EG = c(EG,pm)
	}
	return (EG)
}

#' @describeIn get_EG but \code{get_EG2} computes the Conditional Variance value of the genotype.
get_EG2 <- function( M, p, rdv ){
	L = length(rdv)
	S = 0
	EG = NULL
	g = c(0,1,2)
	for (i in 1:L){
		m = NULL
		for  (j in 1:3){
			L = 1
			for (kk in 1:rdv[i]){
				L = L*M[S + kk,j]
			}
			m = c(m,L*p[j])
		}
		S = S + rdv[i]
		pm = sum(m/sum(m)*g^2) - sum(m/sum(m)*g)^2
		EG = c(EG,pm)
	}	
	return (EG)
}


