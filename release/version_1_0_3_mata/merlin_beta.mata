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

`RM' merlin_logl_beta(`gml' gml)
{
	y 		= merlin_util_depvar(gml)
	linpred = merlin_util_xzb(gml)
	mu 		= invlogit(linpred)
	s 		= asarray(gml.distancb,(gml.model,1))
	mus		= mu:*s
	return(lngamma(s) :- lngamma(mus) :- lngamma(s:-mus) :+ (mus:-1):*log(y) :+ (s:-mus:-1):*log(1:-y))
}

`RM' merlin_beta_expval(`gml' gml, | `RC' t)
{
	if (args()==1) 	return(invlogit(merlin_util_xzb(gml)))
	else 			return(invlogit(merlin_util_xzb(gml,t)))
}

end
