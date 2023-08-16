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

`RM' merlin_logh_hessian_clp(`gml' gml)
{
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
		hazq	= y[,1] :/ 2 :* J(nobs,1,gq[,2]')
		for (q=1;q<=Ngq;q++) {
			X = merlin_util_xz(gml,qp[,q])
			dch = dch :+ exp(merlin_util_xzb(gml,qp[,q])) :* hazq[,q] :* X[,index1] :* X[,index2]
		}
		res = -dch
	}
	else {
		X 	= merlin_util_xz(gml)
		res = - X[,index1] :* X[,index2] :* exp(merlin_util_xzb(gml)) :* y[,1]
	}
	return(res)
}

end
