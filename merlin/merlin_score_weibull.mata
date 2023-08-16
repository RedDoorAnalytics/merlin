*! version 1.0.0 ?????2016

local gml 		struct merlin_struct scalar
local pgml		pointer(struct merlin_struct scalar) scalar
local TR 		transmorphic
local RS 		real scalar
local RC 		real colvector
local SS 		string scalar
local PS 		pointer scalar
local RR 		real rowvector
local RM 		real matrix
local PC 		pointer colvector
local PM 		pointer matrix
local SC 		string colvector

version 14.1

mata:

`RM' merlin_weibull_score_loglambda(`gml' gml, `RS' Xindex)
{
	//core
	gml.survind = 0
	Gx 			= J(gml.N,gml.ndim[gml.Nrelevels],0)	
	y 			= merlin_util_depvar(gml)
	x  			= merlin_util_xz_col_simple(gml,Xindex)
	gam			= merlin_util_dap(gml,1)

	//exactly observed events and/or right censoring -> survival function
	//this first for dimensions
	gml.survind = 2
	Nobs2 		= merlin_get_nobs(gml)
	if (Nobs2) {
		index2 	= merlin_get_surv_index(gml)
		
// 		if (hast | hasbh[2]) {
// 			Ngq 	= gml.chip
// 			chq2 	= J(Nobs2,Ngq,0)
// 			gq 		= merlin_gq(Ngq,"legendre")
// 			if (!haslt) qp2 = y[index2,1] :/ 2 :* J(Nobs2,1,gq[,1]') :+ y[index2,1]:/2
// 			else 		qp2 = (y[index2,1] :- y[index2,3]) :/ 2 :* J(Nobs2,1,gq[,1]') :+ (y[index2,1] :+ y[index2,3]):/2
// 			for (q=1;q<=Ngq;q++) {
// 				chq2[,q] = merlin_util_xzb(gml,qp2[,q]) 
// 			}
// 			if (hasbh[2]) chq2 = chq2 :+ log(merlin_util_bhazard(gml))
// 			chq2 = exp(chq2 :+ log(gam) :+ (gam-1) :* log(qp2))
// 			if (!haslt)	logl[index2] = logl[index2] :- y[index2,1]:/2 :* (chq2 * gq[,2]) 
// 			else 		logl[index2] = logl[index2] :- (y[index2,1]:-y[index2,3]):/2 :* (chq2 * gq[,2]) 
// 		}
// 		else {
			Gx[index2,] = - x[index2] :* exp(merlin_util_xzb(gml)) :* y[index2,1] :^ gam
// 			if (haslt) {
// 				gml.survind 	= 4
// 				index4 			= merlin_get_surv_index(gml)
// 				logl[index4] 	= logl[index4] :+ exp(xb[index4,]) :* y[index4,3] :^ gam
// 			}
// 		}
	}
		
	gml.survind = 1
	if (merlin_get_nobs(gml)) {
		index1 		= merlin_get_surv_index(gml)
		Gx[index1,] = Gx[index1,] :+ x[index1]
	}	
		
	return(Gx)	
		
// 	if (gml.NI[gml.model]) {
// 		nobs 	= gml.Nobs[gml.Nlevels,gml.model]
// 		dch 	= J(nobs,1,0)
// 		Ngq 	= 30
// 		gq 		= merlin_gq(Ngq,"legendre")
// 		qp		= y[,1] :/ 2 :* J(nobs,1,gq[,1]') :+ y[,1]:/2
// 		hazq	= y[,1] :/ 2 :* J(nobs,1,gq[,2]') :* qp :^ (gam-1) :* gam 
// 		for (q=1;q<=Ngq;q++) {
// 			dch = dch :+ exp(merlin_util_xzb(gml,qp[,q])) :* hazq[,q] :* merlin_util_xz(gml,qp[,q])
// 		}
// 		return(y[,2] :* merlin_util_xz(gml) :- dch)
// 	}
// 	else return(merlin_util_xz(gml)[,Xindex] :* (y[,2] :- exp(merlin_util_xzb(gml)) :* y[,1] :^ gam))
}

`RM' merlin_weibull_score_loggamma(`gml' gml, `RS' ind)
{
	gml.survind = 0
	gam	= merlin_util_dap(gml,1)
	y	= merlin_util_depvar(gml)

	if (gml.NI[gml.model]) {
		nobs 	= gml.Nobs[gml.Nlevels,gml.model]
		dch 	= J(nobs,1,0)
		Ngq 	= 30
		gq 		= merlin_gq(Ngq,"legendre")
		qp		= y[,1] :/ 2 :* J(nobs,1,gq[,1]') :+ y[,1]:/2
		hazq	= y[,1]:/2 :* J(nobs,1,gq[,2]') :* qp :^ (gam-1) :* gam :* (gam:*log(qp) :+ 1)
		for (q=1;q<=Ngq;q++) {
			dch = dch :+ exp(merlin_util_xzb(gml,qp[,q])) :* hazq[,q]
		}
		return(y[,2] :* (1 :+ gam :* log(y[,1])) :- dch)
	}
	else return(y[,2] :+ gam :*log(y[,1]) :* (y[,2] :- exp(merlin_util_xzb(gml)) :* y[,1]:^gam))
}



end
