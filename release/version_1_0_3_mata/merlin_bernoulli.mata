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

`RM' merlin_logl_bernoulli(`gml' gml)
{
	y 		= merlin_util_depvar(gml)
	expval 	= (*gml.invlinks[gml.model])(merlin_util_xzb(gml))
	return(y :* log(expval) :+ (1:-y) :* log(1:-expval))
}

`RM' merlin_bernoulli_expval(`gml' gml , | `RC' t)
{
	if (args()==1) 	return((*gml.invlinks[gml.model])(merlin_util_xzb(gml)))
	else 			return((*gml.invlinks[gml.model])(merlin_util_xzb(gml,t)))
}

end
