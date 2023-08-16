*! version 1.0.0 ?????2016

local GMLS 		struct merlin_struct scalar
local pGMLS		pointer(struct merlin_struct scalar) scalar
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

`RM' merlin_p_lognormal_logh(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(merlin_p_lognormal_logpdf(gml,t):-log(merlin_p_lognormal_s(gml,t)))
}

`RM' merlin_p_lognormal_h(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(exp(merlin_p_lognormal_logpdf(gml,t):-log(merlin_p_lognormal_s(gml,t))))
}

`RM' merlin_p_lognormal_ch(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(-log(merlin_p_lognormal_s(gml,t)))
}

`RM' merlin_p_lognormal_logch(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(log(-log(merlin_p_lognormal_s(gml,t))))
}

`RM' merlin_p_lognormal_s(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	mu 		= merlin_util_xzb(gml)
	sigma	= merlin_util_dap(gml,1)

	return(1 :- normal((log(t):-mu):/sigma))
}

`RM' merlin_p_lognormal_pdf(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)

	mod		= gml.model
	mu 		= merlin_util_xzb(gml)
	sigma	= merlin_util_dap(gml,1)

	return(normalden(log(t),mu,sigma) :/ t)
}

`RM' merlin_p_lognormal_logpdf(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)

	mu 		= merlin_util_xzb(gml)
	sigma	= merlin_util_dap(gml,1)
	
	return(lnnormalden(log(t),mu,sigma) :- log(t))
}


end
