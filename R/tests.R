#' Regular Score/Trend test
#'
#' The test_EG function is a regular score/trend test evaluated by asymptotic distribution.
#' 
#' The inputs are two matrices of genotypes for cases and controls.  The matrix should have dimensions: 
#' number of samples (\code{ncase} or \code{ncont}) by number of variants \code{J}.
#'
#' \code{test_EG} - Regular Score/Trend test (asymptotic distbn).\cr
#' \code{test_EGPv} - Regular Score/Trend test (permutation distbn).
#'
#' @param M1 matrix with dimensions number of cases by J variants (double)
#' @param M2 matrix with dimensions number of controls by J variants (double)
#' @return vector of p-values for J variants (double)
test_EG <- function( M1, M2 ){
	L=length(M1[1,])
	cc = rep(0,L)
	for (i in 1:L){
		X = c(M1[,i],M2[,i])
		p = length(M1[,i])/length(X)
		q = 1 - p
		Y = c(rep(1,length(M1[,i])),rep(0,length(M2[,i])))
		Y = Y[!is.na(X)]
		X = X[!is.na(X)]
		a = (glm(Y~1,family='binomial'))
		b = glm(Y~X,family='binomial')
		cc[i] = 1-pchisq(anova(a,b,test='Rao')$Rao[2],1)
	}
	return (cc)
}

#' The test_EGPv function is a regular score/trend test evaluated by permutation distribution.
#' @describeIn test_EG
#' @param P a J by 3 matrix (double) of estimated genotype frequency for J variants: P(G=0), P(G=1) and P(G=2).
#' @param perm number of permutations (int)
test_EGPv <- function( M1, M2, P, perm){
	case = length(M1[,1])
	cont = length(M2[,1])
	L=length(M1[1,])
	cc = rep(0,L)
	for (i in 1:L){
		if (i %% 100 == 1) {cat(i,'\n')}
		X = c(M1[,i],M2[,i]) 
		Y = c(rep(1,length(M1[,i])),rep(0,length(M2[,i])))
		Y = Y[!is.na(X)]
		X = X[!is.na(X)]
		case = sum(Y==1)
		cont = sum(Y==0)
		test = calct(X[1:case],X[(case+1):(case+cont)])
		C = NULL
		n = perm
		for (j in 1:n){
			k = sample(length(X))
			X = X[k]
			a = calct(X[1:case],X[(case+1):(case+cont)])
			C = c(C,a)
		}
		cc[i] = sum(C>=test)/n
	}
	return (cc)
}









#' Robust Variance Score (RVS)
#'
#' The test_EGC function is the RVS test evaluated by asymptotic distribution.
#'
#' RVS uses the Robust Variance Estimates described in the journal Bioinformatics: 
#' Vol.30 no.15 2014, pp.2179 to 2188.
#'
#' The input values are matrices of expected values of genotypes given sequence data for cases and controls.
#'
#' \code{test_EGC} - RVS (asymptotic distbn).\cr
#' \code{test_EGCV} - RVS using Regular Variance Estimates for cases/ctrls (asymptotic distn).\cr
#' \code{test_EGB} - RVS using Robust Variance Estimates for cases/ctrls (bootstrap perm).\cr
#' \code{test_EGBnew} - RVS using the Regular Variance Estimates for cases/ctrls (bootstrap distn and perm). 
#'
#' @param M1 matrix with dimensions of number of cases by J variants (double)
#' @param M2 matrix with dimensions of number of controls by J variants (double)
#' @param P matrix with dimensions J by 3 that consists of estimated genotype frequency for J variants: P(G=0), P(G=1) and P(G=2). (double)
#' @return vector of p-values for J variants (double)
test_EGC <- function( M1, M2, P ){
	L=length(M1[1,])
	cc = rep(0,L)
	for (i in 1:L){
		X = c(M1[,i],M2[,i])
		mx = mean(X)
		p = length(M1[,i])/length(X)
		q = 1 - p
		Y = c(rep(1,length(M1[,i])),rep(0,length(M2[,i])))
		Y = Y[!is.na(X)]
		X = X[!is.na(X)]
		a = (glm(Y~1,family='binomial'))
		b = glm(Y~X,family='binomial')
		p = length(X[Y==1])/length(X)
		q = 1 - p 
		v = q*vp(P[i,]) + p*var(X[Y==0],na.rm=TRUE)
		x = anova(a,b,test='Rao')$Rao[2]*var(X,na.rm=TRUE)/v
		cc[i] = 1-pchisq(x,1)
	}
	return (cc)
}

#' The test_EGCV function is the RVS test evaluated by asymptotic distribution using the Regular Variance Estimates for cases and controls.
#' @describeIn test_EGC 
test_EGCV <- function( M1, M2, P ){
	L=length(M1[1,])
	cc = rep(0,L)
	for (i in 1:L){
		X = c(M1[,i],M2[,i])
		Y = c(rep(1,length(M1[,i])),rep(0,length(M2[,i])))
		Y = Y[!is.na(X)]
		X = X[!is.na(X)]
		a = (glm(Y~1,family='binomial'))
		b = glm(Y~X,family='binomial')
		p = length(X[Y==1])/length(X)
		q = 1 - p 
		v = q*var(X[Y==1]) + p*var(X[Y==0])
		x = anova(a,b,test='Rao')$Rao[2]*var(X)/v
		cc[i] = 1-pchisq(x,1)
	}
	return (cc)
}

#' The test_EGB function is the RVS test evaluated by bootstrap permutation using the Robust Variance Estimates for cases and controls.
#' @describeIn test_EGC
#' @param perm number of permutations (int)
test_EGB <- function( M1, M2, P, perm){
	L=length(M1[1,])
	cc = rep(0,L)
	for (i in 1:L){
		if (i %% 100 == 1) {cat(i,'\n')}
		X = c(M1[,i],M2[,i]) 
		Y = c(rep(1,length(M1[,i])),rep(0,length(M2[,i])))
		Y = Y[!is.na(X)]
		X = X[!is.na(X)]
		case = sum(Y==1)
		cont = sum(Y==0)
		p = length(X[Y==1])/length(X)
		q = 1 - p 
		vcase = vp(P[i,])
		#vcase = var(M1[,i])
		vcont = var(X[Y==0])
		Tobs = (mean(X[Y==1]) - mean(X[Y==0]))/sqrt(vcase/case+vcont/cont) 
		X1 = X[Y==1] - mean(X[Y==1])
		X2 = X[Y==0] - mean(X[Y==0])
		C = NULL
		n = perm
		for (j in 1:n){
			Xca = sample(X1[],case,replace=TRUE)
			Xco = sample(X2[],cont,replace=TRUE)
			vcase = var(Xca)
			vcont = var(Xco)
			C =c(C,(mean(Xca) - mean(Xco))/sqrt(vcase/case+vcont/cont))
		}
		cc[i] = sum(abs(C)>=abs(Tobs))/n
	}
	return (cc)
}

#' The test_EGBnew function is the RVS test evaluated by bootstrap distribution using the Regular Variance Estimates for cases and controls.
#' @describeIn test_EGC
test_EGBnew <- function( M1, M2, P, perm){
	case = length(M1[,1])
	cont = length(M2[,1])
	L=length(M1[1,])
	cc = rep(0,L)
	for (i in 1:L){
		if (i %% 100 == 1) {cat(i,'\n')}
		X = c(M1[,i],M2[,i]) 
		Y = c(rep(1,length(M1[,i])),rep(0,length(M2[,i])))
		Y = Y[!is.na(X)]
		X = X[!is.na(X)]
		p = length(X[Y==1])/length(X)
		q = 1 - p 
		#vcase = vp(P[i,])
		vcase = var(X[Y==1])
		vcont = var(X[Y==0])
		Tobs = (mean(X[Y==1]) - mean(X[Y==0]))/sqrt(vcase/case+vcont/cont) 
		X1 = X[Y==1] - mean(X[Y==1])
		X2 = X[Y==0] - mean(X[Y==0])
		C = NULL
		n = perm
		for (j in 1:n){
			Xca = sample(X1[],case,replace=TRUE)
			Xco = sample(X2[],cont,replace=TRUE)
			vcase = var(Xca)
			vcont = var(Xco)
			C =c(C,(mean(Xca) - mean(Xco))/sqrt(vcase/case+vcont/cont))
		}
		cc[i] = sum(abs(C)>=abs(Tobs))/n
	}
	return (cc)
}


#' Variance of genotypes
#'
#' This function computes the variance of a genotype given the 3 estimated genotype frequencies for variant j: P(G=0), P(G=1) and P(G=2).
#' 
#' vp is called from the \code{test_EGC} and \code{test_EGB} functions.
#'
#' @param P vector of 3 values for variant j, P(G=0), P(G=1) and P(G=2). (double)
#' @return vp variance value (double)
#' @seealso \code{\link{test_EGC}} and \code{\link{test_EGB}}  
vp <- function( P ){
	Sq = 4*P[3] + P[2]
	Sm = 2*P[3] + P[2]
	S = Sq - Sm^2
	return (S)
}


#' Calculate Score Test
#' 
#' This function is called from the \code{test_EGPv} function in the permutation step.
#'
#' @param M1 matrix with dimensions number of cases by J variants (double)
#' @param M2 matrix with dimensions number of controls by J variants (double)
#' @return calct test stat (double)
#' @seealso \code{\link{test_EGPv}} 
calct <- function( M1, M2 ){
	X = c(M1,M2)
	Y = c(rep(1,length(M1)),rep(0,length(M2)))
	p = length(M1)/length(X)
	q = 1 - p
	S = q*sum(M1)-p*sum(M2)
	vs = p*q*(length(X))*var(X)
	x = S^2/vs
	return (x)
}

