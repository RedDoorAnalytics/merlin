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

`RM' merlin_weibull_hessian_ll_ll(`gml' gml)
{
	gam	= merlin_util_dap(gml,1)
	y	= merlin_util_depvar(gml)
	
	nb = asarray(gml.NHbs,gml.model)[1]
	
	index1 = index2 = J(1,0,.)
	for (i=1; i<=nb; i++) {
		refind = 1
		while (refind<=i) {
			index1 = index1,i
			index2 = index2,refind
			refind++
		}
	}
	
	if (gml.NI[gml.model]) {
		nobs 	= gml.Nobs[gml.Nlevels,gml.model]
		dch 	= J(nobs,1,0)
		Ngq 	= 30
		gq 		= merlin_gq(Ngq,"legendre")
		qp		= y[,1] :/ 2 :* J(nobs,1,gq[,1]') :+ y[,1]:/2
		hazq	= y[,1] :/ 2 :* J(nobs,1,gq[,2]') :* qp :^ (gam-1) :* gam 
		for (q=1;q<=Ngq;q++) {
			X = merlin_util_xz(gml,qp[,q])
			dch = dch :+ exp(merlin_util_xzb(gml,qp[,q])) :* hazq[,q] :* X[,index1] :* X[,index2]
		}
		res = -dch
	}
	else {
		X 	= merlin_util_xz(gml)
		res = - X[,index1] :* X[,index2] :* exp(merlin_util_xzb(gml)) :* y[,1] :^ gam
	}
	return(res)
}

`RM' merlin_weibull_hessian_loggamma(`gml' gml)
{
	gam	= merlin_util_dap(gml,1)
	y	= merlin_util_depvar(gml)

	if (gml.NI[gml.model]) {
		nobs 	= gml.Nobs[gml.Nlevels,gml.model]
		dch 	= J(nobs,1,0)
		Ngq 	= 30
		gq 		= merlin_gq(Ngq,"legendre")
		qp		= y[,1] :/ 2 :* J(nobs,1,gq[,1]') :+ y[,1]:/2
		logqp	= log(qp)
		hazq	= y[,1]:/2 :* J(nobs,1,gq[,2]') :* logqp :* gam :* qp:^(gam:-1) :* ((gam:^2 :* logqp :+ 2:*gam) :+ gam  :+ 1:/logqp)
		for (q=1;q<=Ngq;q++) {
			dch = dch :+ hazq[,q] :* exp(merlin_util_xzb(gml,qp[,q])) 
		}
		res = y[,2] :* gam :* log(y[,1]) :- dch
	}
	else {
		logy = log(y[,1])
		res = gam :* logy :* (y[,2] :- exp(merlin_util_xzb(gml)) :* y[,1]:^gam :* (gam :* logy :+ 1) )
	}
	return(res)
}

`RM' merlin_weibull_hessian_ll_lg(`gml' gml)
{
	gam	= merlin_util_dap(gml,1)
	y	= merlin_util_depvar(gml)
	
	if (gml.NI[gml.model]) {
		nobs 	= gml.Nobs[gml.Nlevels,gml.model]
		dch 	= J(nobs,1,0)
		Ngq 	= 30
		gq 		= merlin_gq(Ngq,"legendre")
		qp		= y[,1] :/ 2 :* J(nobs,1,gq[,1]') :+ y[,1]:/2
		hazq	= y[,1] :/ 2 :* J(nobs,1,gq[,2]') :* qp :^ (gam-1) :* gam :* (gam:*log(qp) :+ 1)
		for (q=1;q<=Ngq;q++) {
			dch = dch :+ exp(merlin_util_xzb(gml,qp[,q])) :* hazq[,q] :* merlin_util_xz(gml,qp[,q])
		}
		res = - dch
	}
	else res = - gam :* log(y[,1]) :* merlin_util_xz(gml) :* exp(merlin_util_xzb(gml)) :* y[,1] :^ gam
	return(res)
}

end
