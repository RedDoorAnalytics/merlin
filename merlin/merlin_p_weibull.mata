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

`RM' merlin_p_weibull_logh(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	gam = asarray(gml.distancb,(gml.model,1))
	
	return(merlin_util_xzb(gml,t) :+ log(gam) :+ (gam-1) :* log(t))
}

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
		N = gml.Nobs[gml.Nlevels,gml.model]
		Nq = 30
		gq 	= merlin_gq(Nq,"legendre")
		qp	= t :/ 2 :* J(N,1,gq[,1]') :+ t:/2
		result = J(N,1,0)
		for (q=1; q<=Nq; q++) {
			result = result :+ merlin_p_weibull_h(gml,qp[,q]) :* gq[q,2] :* t:/2
		}
		return(result)
	}
	else return(exp(merlin_util_xzb(gml,t)) :* (t :^ asarray(gml.distancb,(gml.model,1))))
}

`RM' merlin_p_weibull_logch(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)

	if (gml.NI[gml.model]) {	
		N = gml.Nobs[gml.Nlevels,gml.model]
		Nq = 30
		gq 	= merlin_gq(Nq,"legendre")
		qp	= t :/ 2 :* J(N,1,gq[,1]') :+ t:/2
		result = J(N,1,0)
		for (q=1; q<=Nq; q++) {
			result = result :+ merlin_p_weibull_h(gml,qp[,q]) :* gq[q,2] :* t:/2
		}
		return(log(result))
	}
	else return(merlin_util_xzb(gml,t) :+ log(t) :* asarray(gml.distancb,(gml.model,1)))
}

`RM' merlin_p_weibull_s(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	return(exp(-merlin_p_weibull_ch(gml,t)))
}

`RM' merlin_p_weibull_f(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	return(1:-merlin_p_weibull_s(gml,t))
}

`RM' merlin_p_weibull_rmst(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	N = gml.Nobs[gml.Nlevels,gml.model]
	Nq = 30
	gq 	= merlin_gq(Nq,"legendre")
	qp	= t :/ 2 :* J(N,1,gq[,1]') :+ t:/2

	result = J(N,1,0)
	for (q=1; q<=Nq; q++) {
		result = result :+ merlin_p_weibull_s(gml,qp[,q]) :* gq[q,2] :* t:/2
	}
	return(result)
}

`RM' merlin_p_weibull_rmft(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(t:-merlin_p_weibull_rmst(gml,t))
}

`RM' merlin_p_weibull_mu(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	gam = asarray(gml.distancb,(gml.model,1))
	return(gamma(1+1/gam) :* exp(-merlin_util_xzb(gml,t):/gam))	
}

`RM' merlin_p_weibull_dens(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(merlin_p_weibull_h(gml,t):*merlin_p_weibull_s(gml,t))
}

end
