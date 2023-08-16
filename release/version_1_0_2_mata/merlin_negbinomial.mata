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

`RM' merlin_logl_negbinomial(`gml' gml)
{
	`RS' m1, alpha1
	`RM' p
	`RC' y
	y = merlin_util_depvar(gml)	
	linpred = merlin_util_xzb(gml)
	alpha1 = exp(asarray(gml.distancb,(gml.model,1)))
	m1 = 1/alpha1
	p = 1:/(1:+ alpha1:*log(linpred))
	return(lngamma(y:+m1) :- lngamma(y:+1) :- lngamma(m1) :+ m1 :* log(p) :+ y:*log(1:-p))
}

`RM' merlin_negbinomial_expval(`gml' gml, | `RC' t)
{
	if (args()==1) 	return(exp(merlin_util_xzb(gml)))
	else 			return(exp(merlin_util_xzb(gml,t)))
}

end
