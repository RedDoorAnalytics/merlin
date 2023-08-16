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

`RM' merlin_logch_score_clp(`gml' gml)
{
	y	= merlin_util_depvar(gml)	
	return(merlin_util_xz(gml,y[,1]) :* (y[,2] :- exp(merlin_util_xzb(gml,y[,1]))) :+ y[,2] :* merlin_util_xz_deriv(gml,y[,1]) :/ merlin_util_xzb_deriv(gml,y[,1]))
}

end
