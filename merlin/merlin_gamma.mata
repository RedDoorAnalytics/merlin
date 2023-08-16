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

`RM' merlin_logl_gamma(`gml' gml , | `RM' G, `RM' H)
{
	y 		= merlin_util_depvar(gml)
	mu 		= exp(merlin_util_xzb(gml))
	v 		= 1:/asarray(gml.distancb,(gml.model,1)):^2	
	logl 	= (v:-1):* log(y) :- v:/mu:*y :+ v :* log(v) :- v :* log(mu) :- lngamma(v)
	return(logl)
}

`RM' merlin_gamma_expval(`gml' gml, | `RC' t)
{
	if (args()==1)	return(exp(merlin_util_xzb(gml))) 
	else 			return(exp(merlin_util_xzb(gml,t))) 
}

end

