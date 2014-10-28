#' EM calculation
#'
#' Given a n by 3 matrix, M containing likelihoods this function uses the EM algorithm to estimate genotype frequencies.
#'
#' This part is the new part of the analysis where we can use reads instead of hard calls.\cr
#' For this to work we need to estimate genotype frequencies.
#'
#' Data consists of matrix n by 3\cr
#' p = P(G=0)\cr
#' q = P(G=1)\cr
#'
#' @param M matrix of likelihoods derived from PL data.
#' @return a vector of three values: p0, q0, 1-p0-q0
EM_calc <- function(M){
	p_0 = 0.15
	q_0 = 0.15
	q_n = 1
	p_n = 0
	d_n = 0
	k = 0
	while ((p_n - p_0)^2 + (q_n - q_0)^2>0.000001){
		d_0 = 1-p_0 - q_0
		v = c(p_0,q_0,d_0)
		p_D = M%*%(v)
		E_p = M[,1]*p_0/p_D
		E_q = M[,2]*q_0/p_D
		p_n= p_0
		q_n = q_0
		d_n = 1-q_0-p_0
		p_0 = sum(E_p)/length(E_p)
		q_0 = sum(E_q)/length(E_q)
		k = k+1
		if (k==1000){
			cat('hi','\n')
			return ( c(p_0, q_0, 1-p_0-q_0) )
		}
	}	
	return ( c(p_0, q_0, 1-p_0-q_0) )
}


