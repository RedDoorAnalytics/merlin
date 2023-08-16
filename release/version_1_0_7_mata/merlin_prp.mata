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

`RM' merlin_logl_prp(`gml' gml)
{
	mod		= gml.model
	y 		= merlin_util_depvar(gml)
	brcs	= asarray(gml.distancb,(mod,1))
	logch	= merlin_util_xzb(gml) :+ asarray(gml.distancb,(mod,2)) * brcs
	drcs 	= asarray(gml.distancb,(mod,6)) * brcs

	//penalty
	
		knots 		= asarray(gml.distancb,(gml.model,3))
		hasorthog 	= asarray(gml.distancb,(gml.model,4))
		Ngq			= 100
		gq 			= merlin_gq(Ngq,"legendre")
		mt			= minmax(log(y))
		nodes 		= (mt[2] :- mt[1]) :* gq[,1]:/2 :+ (mt[2] :+ mt[1]):/2
		weights 	= (mt[2] :- mt[1]) :* gq[,2]':/2
		
		lambda = 1
		
		if (hasorthog) {
			rmat 	= asarray(gml.distancb,(gml.model,5))
			pen 	= 0
			sp		= merlin_rcs(nodes,knots,2,rmat) //* brcs
			for (k=1; k<=Ngq; k++) {
				pen = pen :+ weights[k] :* brcs' * (sp[k,]'*sp[k,]) * brcs
			} 				
		}
		else {
			pen 	= J(rows(y),1,0)
			for (k=1; k<=Ngq; k++) {
				sp	= merlin_rcs(log(nodes[k]),knots,2) * brcs 
				pen = pen :+ weights[k] :* sp:^2
			} 
		}
	
		pen = -pen :* lambda :/2 :/rows(y)
		
	//logl
	if (gml.nobhaz[mod]) {
		return(pen :+ y[,2] :* (logch :+ log(drcs :+ y[,1] :* merlin_util_xzb_deriv(gml,y[,1]))) :- exp(logch))
	}
	else {
		return(pen :+ y[,2] :* log(exp(logch :+ log(drcs:/y[,1] :+ merlin_util_xzb_deriv(gml,y[,1]))) :+ asarray(gml.bhazards,mod)) :- exp(logch))
	}
}

`RM' merlin_prp_logch(`gml' gml)
{
	return(merlin_util_xzb(gml) :+ asarray(gml.distancb,(gml.model,2)) * asarray(gml.distancb,(mod,1)))
}

`RM' merlin_prp_s(`gml' gml)
{
	return(exp(-exp(merlin_prp_logch(gml))))
}

`RM' merlin_prp_logh(`gml' gml)
{
	mod		= gml.model
	y 		= merlin_util_depvar(gml)
	brcs	= asarray(gml.distancb,(mod,1))

	logch	= merlin_util_xzb(gml) :+ asarray(gml.distancb,(mod,2)) * asarray(gml.distancb,(mod,1))
	drcs 	= asarray(gml.distancb,(mod,6)) * brcs

	return(logch :+ log(drcs :+ merlin_util_xzb_deriv(gml,y[,1])):/y[,1])
}

end
