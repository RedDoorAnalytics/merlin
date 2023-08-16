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

/*
-> functions for family(rp, ...)
*/

mata:

`RM' merlin_p_rp_logch(`gml' gml,| `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	brcs 	= asarray(gml.distancb,(gml.model,1))
	knots 	= asarray(gml.distancb,(gml.model,3))
	hasorth = asarray(gml.distancb,(gml.model,4))
	if (hasorth) {
		rmat = asarray(gml.distancb,(gml.model,5))
		return(merlin_util_xzb(gml,t) :+ merlin_rcs(log(t),knots,0,rmat) * brcs )
	}
	else {
		return(merlin_util_xzb(gml,t) :+ merlin_rcs(log(t),knots,0) * brcs )
	}
}

`RM' merlin_p_rp_ch(`gml' gml,| `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(exp(merlin_p_rp_logch(gml,t)))
}

`RM' merlin_p_rp_s(`gml' gml,| `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(exp(-merlin_p_rp_ch(gml,t)))
}

`RM' merlin_p_rp_logh(`gml' gml,| `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	hast	= gml.istimedep[gml.model,1]
	brcs 	= asarray(gml.distancb,(gml.model,1))
	knots 	= asarray(gml.distancb,(gml.model,3))
	hasorth = asarray(gml.distancb,(gml.model,4))
	logt    = log(t)
	if (hasorth) {
		rmat = asarray(gml.distancb,(gml.model,5))
		logh = merlin_util_xzb(gml,t) :+ merlin_rcs(logt,knots,0,rmat) * brcs 
		if (hast) 	logh = logh :+ log((merlin_rcs(logt,knots,1,rmat) * brcs) :/ t :+ merlin_util_xzb_deriv(gml,t))
		else 		logh = logh :+ log((merlin_rcs(logt,knots,1,rmat) * brcs) :/ t)
	}
	else {
		logh = merlin_util_xzb(gml,t) :+ merlin_rcs(logt,knots,0) * brcs 
		if (hast) 	logh = logh :+ log((merlin_rcs(logt,knots,1) * brcs) :/ t:+ merlin_util_xzb_deriv(gml,t))
		else 		logh = logh :+ log((merlin_rcs(logt,knots,1) * brcs) :/ t)
	}
	return(logh)
}

`RM' merlin_p_rp_h(`gml' gml,| `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(exp(merlin_p_rp_logh(gml,t)))
}


`RM' merlin_p_rp_f(`gml' gml,| `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(1 :- merlin_p_rp_s(gml,t))
}

`RM' merlin_p_rp_dens(`gml' gml,| `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(merlin_p_rp_h(gml,t):*merlin_p_rp_s(gml,t))
}

end
