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

`RM' merlin_p_ggamma_logh(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(merlin_p_ggamma_logpdf(gml,t):-log(merlin_p_ggamma_s(gml,t)))
}

`RM' merlin_p_ggamma_h(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(exp(merlin_p_ggamma_logpdf(gml,t):-log(merlin_p_ggamma_s(gml,t))))
}

`RM' merlin_p_ggamma_ch(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(-log(merlin_p_ggamma_s(gml,t)))
}

`RM' merlin_p_ggamma_logch(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(log(-log(merlin_p_ggamma_s(gml,t))))
}

`RM' merlin_p_ggamma_s(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	mod		= gml.model
	mu 		= merlin_util_xzb(gml)
	sigma	= merlin_util_dap(gml,1)
	kapp	= merlin_util_dap(gml,2)
	gamm	= abs(kapp)^(-2)
	
	z		= sign(kapp) :* (log(t) :- mu):/sigma
	u		= gamm :* exp(abs(kapp) :* z)
	
	if (kapp>0) 		return(1 :- gammap(gamm,u))
	else if (kapp<0) 	return(gammap(gamm,u))
	else 				return(1 :- normal(z))
}

`RM' merlin_p_ggamma_pdf(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)

	mod	= gml.model
	mu 	= merlin_util_xzb(gml)
	sigma	= merlin_util_dap(gml,1)
	kapp	= merlin_util_dap(gml,2)
	gamm	= abs(kapp)^(-2)

	z	= sign(kapp) :* (log(t) :- mu):/sigma
	u	= gamm :* exp(abs(kapp) :* z)

	if (kapp==0) 	return(exp(-log(sigma :* t :* sqrt(2*pi())) :- (z:^2):/2))
	else 		return(exp(gamm:*log(gamm) :- log(sigma :* t :* sqrt(gamm) :* gamma(gamm)) :+ z:*sqrt(gamm) :- u))	
}

`RM' merlin_p_ggamma_logpdf(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)

	mod	= gml.model
	mu 	= merlin_util_xzb(gml)
	sigma	= merlin_util_dap(gml,1)
	kapp	= merlin_util_dap(gml,2)
	gamm	= abs(kapp)^(-2)

	z	= sign(kapp) :* (log(t) :- mu):/sigma
	u	= gamm :* exp(abs(kapp) :* z)

	if (kapp==0) 	return(-log(sigma :* t :* sqrt(2*pi())) :- (z:^2):/2)
	else 		return(gamm:*log(gamm) :- log(sigma :* t :* sqrt(gamm) :* gamma(gamm)) :+ z:*sqrt(gamm) :- u)
}


end
