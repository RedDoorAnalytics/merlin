*! version 1.0.0 ?????2016

local GMLS 		struct merlin_struct scalar
local pGMLS		pointer(struct merlin_struct scalar) scalar
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

`RM' merlin_p_weibull_h(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	gam = asarray(gml.distancb,(gml.model,1))
	
	return(exp(merlin_util_xzb(gml,t)) :* gam :* t :^ (gam - 1))
}

`RM' merlin_p_weibull_ch(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	if (gml.NI[gml.model]) {	
		gq 	= merlin_gq(15,"legendre")
		qp	= t :/ 2 :* J(gml.N,1,gq[,1]') :+ t:/2
		result = J(gml.N,1,0)
		for (q=1; q<=15; q++) {
			result = result :+ merlin_p_weibull_h(gml,qp[,q]) :* gq[q,2] :* t:/2
		}
		return(result)
	}
	else return(exp(merlin_util_xzb(gml,t)) :* (t :^ asarray(gml.distancb,(gml.model,1))))
}

`RM' merlin_p_weibull_s(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	return(exp(-merlin_p_weibull_ch(gml,t)))
}

`RM' merlin_p_weibull_rmst(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	gq 	= merlin_gq(15,"legendre")
	qp	= t :/ 2 :* J(gml.N,1,gq[,1]') :+ t:/2

	result = J(gml.N,1,0)
	for (q=1; q<=15; q++) {
		result = result :+ merlin_p_weibull_s(gml,qp[,q]) :* gq[q,2] :* t:/2
	}
	return(result)
}

`RM' merlin_p_weibull_mu(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	gam = asarray(gml.distancb,(gml.model,1))
	return(gamma(1+1/gam) :* exp(-merlin_util_xzb(gml,t):/gam))	
}

end
