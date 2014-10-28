#' Generate sequence data under null or alt hypothesis
#'
#' This function is designed to mimic low and high read depth data given parameters that simulate hard calls.
#' 
#' @param L number of variants (integer) >0
#' @param ncase number of cases (integer) >0
#' @param ncont number of controls (integer) >0
#' @param mdcase average read depth in cases (double) >0
#' @param sddcase standard deviation of read depth in cases  (double) >0
#' @param mdcont average read depth in controls (double) >0
#' @param sddcont standard deviation of read depth in controls (double) >0  
#' @param em average error rate, probability that sequence call is wrong (double) >0
#' @param esd standard deviation of error rate (double) >0
#' @param mmaf MAF for L variants, HWE is assumed >0
#' @return For null, outputs an object:
#'
#' MM - conditional expected value E(Gij|Dij)\cr
#' P  - estimated genotype frequencies P(G=0), P(G=1), P(G=2)\cr
#' G  - true genotype from which sequence data generated
generate_seqdata_null <- function( L, ncase, ncont, mdcase, sddcase, mdcont, sddcont, em, esd, mmaf ) {
	MAF = NULL
	A1 = 'C'
	A2 = 'T'
	MM = NULL
	G = NULL
	P = NULL
	for (i in 1:L){
		if (i %% 100 == 1) {cat(i,'\n')}
		v = NULL
		erv = NULL
		rdv = NULL
		gv = NULL
		maf = mmaf[i]
		MAF = c(MAF,maf)
		for (j in 1:ncase){
			rd = round(sddcase*rnorm(1) + mdcase)
			if (rd <= 0){rd=1}
			error = esd*rnorm(rd) + em
			k = rbinom(1,2,maf)
			gv = c(gv,k)
			gen_vect = c('CC','CT','TT')
			if (k==2){genotype=c('C','C')}
			if (k==1){genotype=c('C','T')}
			if (k==0){genotype=c('T','T')}
			a = seq_call(genotype, error, gen_vect)
			v = c(v,a)
			erv = c(erv,error)
			rdv = c(rdv,rd)
		}
		for (j in 1:ncont){
			rd = round(sddcont*rnorm(1) + mdcont)
			if (rd <= 0){rd=1}
			error = esd*rnorm(rd) + em
			k = rbinom(1,2,maf)
			gv = c(gv,k)
			gen_vect = c('CC','CT','TT')
			if (k==2){genotype=c('C','C')}
			if (k==1){genotype=c('C','T')}
			if (k==0){genotype=c('T','T')}
			a = seq_call(genotype, error, gen_vect)
			v = c(v,a)
			erv = c(erv,error)
			rdv = c(rdv,rd)
		}
		G = cbind(G,gv)
		M = get_Mr(erv,v,rdv)
		Mm = get_M(erv,v)
		p = EM_calc(M)
		#p = c(0.9^2,2*0.9*0.1,0.1^2)
		P = rbind(P,p)
		if ((sum(p< (-0.00001))>0) || (sum(p>1.00001)<0)){
			cat('Problem!!!')
			return (NULL)
		}
		EG = get_EG(Mm,p,rdv)
		MM = cbind(MM,EG)
	}
	return (list(MM=MM,P=P,G=G))
}


#' @describeIn generate_seqdata_null With MAF for cases/controls.
#' @param mmafCa MAF for L variants in cases, HWE is assumed >0
#' @param mmafCo MAF for L variants in controls, HWE is assumed >0
#' @return For alt, outputs same object with additional:
#'
#' t=MM1 - conditional expected value E(Gij|Dij)
#'
#' Note to Andriy: The seq_call() function that is called from  generate_seqdata_alt(), doesn't have input args: gen_freq and R.
#' I commented these two args out for now.  Also, gen_vect is an empty arg; see the seq_call() code.  --Ted.
#generate_seqdata_alt <- function ( L, ncase, ncont, mdcase, sddcase, mdcont, sddcont, em, esd, mmafCa, mmafCo, R) {
generate_seqdata_alt <- function ( L, ncase, ncont, mdcase, sddcase, mdcont, sddcont, em, esd, mmafCa, mmafCo) {
	MAF = NULL
	A1 = 'C'
	A2 = 'T'
	MM = NULL
	MM1 = NULL
	G = NULL
	P = NULL
	for (i in 1:L){
		if (i %% 100 == 1) {cat(i,'\n')}
		v = NULL
		erv = NULL
		rdv = NULL
		gv = NULL
		maf = mmafCa[i]
		MAF = c(MAF,maf)
		for (j in 1:ncase){
			rd = round(sddcase*rnorm(1) + mdcase)
			if (rd <= 0){rd=1}
			error = esd*rnorm(rd) + em
			k = rbinom(1,2,maf)
			gv = c(gv,k)
			gen_vect = c('CC','CT','TT')
			if (k==2){genotype=c('C','C')}
			if (k==1){genotype=c('C','T')}
			if (k==0){genotype=c('T','T')}
			# a = seq_call(genotype, error, gen_vect, gen_freq,R)
			a = seq_call(genotype, error, gen_vect)
			v = c(v,a)
			erv = c(erv,error)
			rdv = c(rdv,rd)
		}
		maf = mmafCo[i]
		for (j in 1:ncont){
			rd = round(sddcont*rnorm(1) + mdcont)
			if (rd <= 0){rd=1}
			error = esd*rnorm(rd) + em
			k = rbinom(1,2,maf)
			gv = c(gv,k)
			gen_vect = c('CC','CT','TT')
			if (k==2){genotype=c('C','C')}
			if (k==1){genotype=c('C','T')}
			if (k==0){genotype=c('T','T')}
			# a = seq_call(genotype, error, gen_vect, gen_freq,R)
			a = seq_call(genotype, error, gen_vect)
			v = c(v,a)
			erv = c(erv,error)
			rdv = c(rdv,rd)
		}
		G = cbind(G,gv)
		M = get_Mr(erv,v,rdv)
		Mm = get_M(erv,v)
		p = EM_calc(M)
		#p = c(0.9^2,2*0.9*0.1,0.1^2)
		P = rbind(P,p)
		if ((sum(p< (-0.00001))>0) || (sum(p>1.00001)<0)){
			cat('Problem!!!')
			return (NULL)
		}
		EG = get_EG(Mm,p,rdv)
		EG2 = get_EG2(Mm,p,rdv)
		MM = cbind(MM,EG)
		MM1 = cbind(MM1,EG2)
	}
	return (list(MM=MM,P=P,G=G,t=MM1))
}

