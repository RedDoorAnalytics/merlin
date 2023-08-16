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

`RM' merlin_rp_score_clp(`gml' gml)
{
	y	= merlin_util_depvar(gml)
	
	brcs 	= asarray(gml.distancb,(gml.model,1))
	knots 	= asarray(gml.distancb,(gml.model,3))
	hasorth = asarray(gml.distancb,(gml.model,4))
	logt    = log(y[,1])
	if (hasorth) {
		rmat = asarray(gml.distancb,(gml.model,5))
		logh = merlin_util_xz(gml,y[,1]) :+ merlin_util_xz_deriv(gml,y[,1]) :/ ((merlin_rcs(logt,knots,1,rmat) * brcs) :/ y[,1] :+ merlin_util_xzb_deriv(gml,y[,1]))
	}
	else {
		logh = merlin_util_xz(gml,y[,1]) :+ merlin_util_xz_deriv(gml,y[,1]) :/ ((merlin_rcs(logt,knots,1) * brcs) :/ y[,1] :+ merlin_util_xzb_deriv(gml,y[,1]) )
	}
	
	return(y[,2] :* logh :-  merlin_rp_ch(gml,y[,1]) :* merlin_util_xz(gml,y[,1]))
	
}

`RM' merlin_rp_score_dap(`gml' gml)
{
	y	= merlin_util_depvar(gml)
	
	brcs 	= asarray(gml.distancb,(gml.model,1))
	knots 	= asarray(gml.distancb,(gml.model,3))
	hasorth = asarray(gml.distancb,(gml.model,4))
	logt    = log(y[,1])
	if (hasorth) {
		rmat = asarray(gml.distancb,(gml.model,5))
		logh = merlin_rcs(logt,knots,0,rmat) :+ (merlin_rcs(logt,knots,1,rmat) :/ y[,1]) :/ ((merlin_rcs(logt,knots,1,rmat) * brcs) :/ y[,1] :+ merlin_util_xzb_deriv(gml,y[,1]))
		return(y[,2] :* logh :- merlin_rp_ch(gml,y[,1]) :* merlin_rcs(logt,knots,0,rmat))
	}
	else {
		logh = merlin_rcs(logt,knots,0) :+ (merlin_rcs(logt,knots,1) :/ y[,1]) :/ ((merlin_rcs(logt,knots,1) * brcs) :/ y[,1] :+ merlin_util_xzb_deriv(gml,y[,1]) )
		return(y[,2] :* logh :- merlin_rp_ch(gml,y[,1]) :* merlin_rcs(logt,knots,0))
	}	
}


end
