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

`RM' merlin_p_rp_h(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(exp(merlin_p_rp_logh(gml,t)))
}


`RM' merlin_p_rp_logh(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	mod 		= gml.model
	brp 		= asarray(gml.distancb,(mod,1))
	knots 		= asarray(gml.distancb,(mod,3))
	hasorthog 	= asarray(gml.distancb,(mod,4))
	logt		= log(t)
	
	if (hasorthog) {
		rmat 	= asarray(gml.distancb,(mod,5))
		logch 	= merlin_util_xzb(gml,t) :+ merlin_rcs(logt,knots,0,rmat) * brp
		drcs	= merlin_rcs(logt,knots,1,rmat) * brp
	}
	else {
		logch 	= merlin_util_xzb(gml,t) :+ merlin_rcs(logt,knots,0) * brp
		drcs	= merlin_rcs(logt,knots,1) * brp
	}

	return(logch :+ log(drcs:/t :+ merlin_util_xzb_deriv(gml,t))) 
}

`RM' merlin_p_rp_ch(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	mod 		= gml.model
	brp 		= asarray(gml.distancb,(mod,1))
	knots 		= asarray(gml.distancb,(mod,3))
	hasorthog 	= asarray(gml.distancb,(mod,4))
	logt		= log(t)

	if (hasorthog) {
		rmat 	= asarray(gml.distancb,(mod,5))
		return(exp(merlin_util_xzb(gml,t) :+ merlin_rcs(logt,knots,0,rmat) * brp))
	}
	else {
		return(exp(merlin_util_xzb(gml,t) :+ merlin_rcs(logt,knots,0) * brp))
	}

}

`RM' merlin_p_rp_logch(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	mod 		= gml.model
	brp 		= asarray(gml.distancb,(mod,1))
	knots 		= asarray(gml.distancb,(mod,3))
	hasorthog 	= asarray(gml.distancb,(mod,4))
	logt		= log(t)

	if (hasorthog) {
		rmat 	= asarray(gml.distancb,(mod,5))
		return(merlin_util_xzb(gml,t) :+ merlin_rcs(logt,knots,0,rmat) * brp)
	}
	else {
		return(merlin_util_xzb(gml,t) :+ merlin_rcs(logt,knots,0) * brp)
	}
}

`RM' merlin_p_rp_s(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	return(exp(-merlin_p_rp_ch(gml,t)))
}

`RM' merlin_p_rp_f(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	return(1:-merlin_p_rp_s(gml,t))
}

`RM' merlin_p_rp_rmst(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	N = gml.Nobs[gml.Nlevels,gml.model]
	Nq = 30
	gq 	= merlin_gq(Nq,"legendre")
	qp	= t :/ 2 :* J(N,1,gq[,1]') :+ t:/2

	result = J(N,1,0)
	for (q=1; q<=Nq; q++) {
		result = result :+ merlin_p_rp_s(gml,qp[,q]) :* gq[q,2] :* t:/2
	}
	return(result)
}

`RM' merlin_p_rp_rmft(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(t:-merlin_p_rp_rmst(gml,t))
}

end
