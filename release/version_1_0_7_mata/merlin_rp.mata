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
-> functions for family(rp, ...)
*/

`RM' merlin_logl_rp(`gml' gml)
{
	mod		= gml.model
	y 		= merlin_util_depvar(gml)
	brcs	= asarray(gml.distancb,(mod,1))
	logch	= merlin_util_xzb(gml) :+ asarray(gml.distancb,(mod,2)) * brcs
	drcs 	= asarray(gml.distancb,(mod,6)) * brcs

	if (gml.nobhaz[mod]) {
		return(y[,2] :* (logch :+ log(drcs :+ y[,1] :*  merlin_util_xzb_deriv(gml,y[,1]))) :- exp(logch))
	}
	else {
		return(y[,2] :* log(exp(logch :+ log(drcs:/y[,1] :+ merlin_util_xzb_deriv(gml,y[,1]))) :+ asarray(gml.bhazards,mod)) :- exp(logch))
	}
}

`RM' merlin_rp_logch(`gml' gml)
{
	return(merlin_util_xzb(gml) :+ asarray(gml.distancb,(gml.model,2)) * asarray(gml.distancb,(mod,1)))
}

`RM' merlin_rp_s(`gml' gml)
{
	return(exp(-exp(merlin_rp_logch(gml))))
}

`RM' merlin_rp_logh(`gml' gml)
{
	mod		= gml.model
	y 		= merlin_util_depvar(gml)
	brcs	= asarray(gml.distancb,(mod,1))

	logch	= merlin_util_xzb(gml) :+ asarray(gml.distancb,(mod,2)) * asarray(gml.distancb,(mod,1))
	drcs 	= asarray(gml.distancb,(mod,6)) * brcs

	return(logch :+ log(drcs :+ merlin_util_xzb_deriv(gml,y[,1])):/y[,1])
}

end
