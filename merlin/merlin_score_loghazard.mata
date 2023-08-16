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

`RM' merlin_logh_score_clp(`gml' gml, | `RS' Xindex)
{
	y	= merlin_util_depvar(gml)

	if (gml.NI[gml.model]) {
		nobs 	= gml.Nobs[gml.Nlevels,gml.model]
		dch 	= J(nobs,1,0)
		Ngq 	= 30
		gq 		= merlin_gq(Ngq,"legendre")
		qp		= y[,1] :/ 2 :* J(nobs,1,gq[,1]') :+ y[,1]:/2
		hazq	= y[,1] :/ 2 :* J(nobs,1,gq[,2]') 
		for (q=1;q<=Ngq;q++) {
			dch = dch :+ exp(merlin_util_xzb(gml,qp[,q])) :* hazq[,q] :* merlin_util_xz(gml,qp[,q])
		}
		return(y[,2] :* merlin_util_xz(gml) :- dch)
	}
	else {
		if (args()==1) {
			return(merlin_util_xz(gml) :* (y[,2] :- exp(merlin_util_xzb(gml)) :* y[,1]))
		}
		else {
			return(merlin_util_xz(gml)[,Xindex] :* (y[,2] :- exp(merlin_util_xzb(gml)) :* y[,1]))
		}
	}
}

end
