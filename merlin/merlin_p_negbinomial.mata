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

`RM' merlin_p_negbinomial_mu(`gml' gml , | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(exp(merlin_util_xzb(gml,t)))
}

end
