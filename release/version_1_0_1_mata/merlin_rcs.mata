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

/*
-> functions for family(rcs, ...)
*/

`RM' merlin_logl_rcs(`gml' gml)
{
	mod			= gml.model
	y 			= merlin_util_depvar(gml)
	logy		= log(y[,1])
	brcs		= asarray(gml.distancb,(mod,1))
	knots 		= asarray(gml.distancb,(mod,3))
	hasorthog 	= asarray(gml.distancb,(mod,4))
	if (hasorthog) rmat = asarray(gml.distancb,(mod,5))

	logh		= merlin_util_xzb(gml) :+ asarray(gml.distancb,(mod,2)) * brcs
	
	ch 			= J(gml.Nobs[gml.Nlevels,mod],1,0)
	hq			= asarray(gml.haznodes,mod)
	loghq 		= log(hq)
	loghw		= log(asarray(gml.hazweights,mod))

	for (k=1;k<=gml.hazNnodes[mod];k++) {
		loghk = merlin_util_xzb(gml,hq[,k])
		if (hasorthog) 	loghk = loghk :+ merlin_rcs(loghq[,k],knots,0,rmat) 	* brcs
		else 			loghk = loghk :+ merlin_rcs(loghq[,k],knots,0) 		* brcs
		ch = ch :+ exp(loghw[,k] :+ loghk)
	}
	
	if (gml.nobhaz[mod]) return(y[,2] :* (logh :+ logy) :- ch)
	else return(y[,2] :* log(exp(logh) :+ asarray(gml.bhazards,mod)) :- ch)
}

`RM' merlin_rcs_logh(`gml' gml)
{
	return(merlin_util_xzb(gml) :+ asarray(gml.distancb,(gml.model,2)))
}

`RM' merlin_rcs_ch_ni(`gml' gml)
{
	mod			= gml.model
	brcs		= asarray(gml.distancb,(mod,1))
	knots 		= asarray(gml.distancb,(mod,3))
	hasorthog 	= asarray(gml.distancb,(mod,4))
	if (hasorthog) rmat = asarray(gml.distancb,(mod,5))

	ch 			= J(gml.Nobs[gml.Nlevels,mod],1,0)
	hq			= asarray(gml.haznodes,mod)
	loghq 		= log(hq)
	loghw		= log(asarray(gml.hazweights,mod))

	for (k=1;k<=gml.hazNnodes[mod];k++) {
		loghk = merlin_util_xzb(gml,hq[,k])
		if (hasorthog) 	loghk = loghk :+ merlin_rcs(loghq[,k],knots,0,rmat) 	* brcs
		else 			loghk = loghk :+ merlin_rcs(loghq[,k],knots,0) 		* brcs
		ch = ch :+ exp(loghw[,k] :+ loghk)
	}

	return(ch)
}

`RM' merlin_rcs_s_ni(`gml' gml)
{
	return(exp(-merlin_rcs_ch_ni(gml)))
}

end
