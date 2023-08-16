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

`RM' merlin_logl_qtile(`gml' gml)
{
	y 			= merlin_util_depvar(gml)
	mu 			= merlin_util_xzb(gml)
	sigma 		= asarray(gml.distancb,(gml.model,1))
	logsigma 	= log(sigma)
	tau 		= gml.qtiles[gml.model]	
	omtau		= 1-tau
	v 			= y :- mu
	ncol 		= cols(mu)
	row 		= J(rows(mu),ncol,0)

	for (c=1; c<=ncol; c++) {
		index1 			= selectindex(v[,c]:>0)
		row[index1,c] 	= v[index1,c] :* tau
		index2 			= selectindex((-v[,c]):>0) 
		row[index2,c] 	= -v[index2,c] :* omtau
	}

	return(log(tau) :+ log(omtau) :- logsigma :- row :/ sigma)
}


end
