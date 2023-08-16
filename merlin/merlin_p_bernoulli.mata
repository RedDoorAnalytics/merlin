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

`RM' merlin_p_bernoulli_mu(`gml' gml , | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return((*gml.invlinks[gml.model])(merlin_util_xzb(gml,t)))
}

`RM' merlin_p_over_bernoulli_mu(`gml' gml , `PS' pf, | `RC' t)
{
	not = args()==2
	if (not) t = merlin_util_timevar(gml)
	
	mod = gml.model
	overmod = strtoreal(st_local("overmodel"))
	overatvars = tokens(st_local("overatvars"))
	overatvals = strtoreal(tokens(st_local("overatvals")))
	
	//one
		
		//outcome
		//one
		merlin_util_xzb_update(gml,st_global(" e(response"+strofreal(overmod)+")"), 1)
		out1 = (*pf)(gml,t)
		//zero
		merlin_util_xzb_update(gml,st_global(" e(response"+strofreal(overmod)+")"), 0)
		out0 = (*pf)(gml,t)
		
		//overoutcome
		gml.model = overmod
		merlin_util_xzb_update(gml,overatvars,overatvals,overmod)
		prob1 = merlin_p_bernoulli_mu(gml, t)
		gml.model = mod
		
	return(out1 :* prob1 :+ out0 :* (1:-prob1))
	
}

end
