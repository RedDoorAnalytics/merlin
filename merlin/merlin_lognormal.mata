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

`RM' merlin_logl_lognormal(`gml' gml, `RM' G, `RM' H)
{
	gml.survind = 0

	mod		= gml.model
	y 		= merlin_util_depvar(gml)
	mu 		= merlin_util_xzb(gml)
	sigma	= merlin_util_dap(gml,1)

	logl 	= J(merlin_get_nobs(gml,mod),1,0)
	
	//exact events
	gml.survind = 1
	if (merlin_get_nobs(gml)) {
		index 			= merlin_get_surv_index(gml)
		logl[index,] 	= lnnormalden(log(y[index,1]),mu[index,],sigma) :- log(y[index,1])
	}

	//right censored
	gml.survind = 6
	if (merlin_get_nobs(gml)) {
		index 			= merlin_get_surv_index(gml)
		logl[index,] 	= log(1 :- normal((log(y[index,1]):-mu[index,]):/sigma))
	}

	return(logl)
}

`RM' merlin_logl_lognormal_ml(`gml' gml)
{
	gml.survind = 0

	mod		= gml.model
	y 		= merlin_util_depvar(gml)
	mu 		= merlin_util_xzb(gml)
	sigma	= merlin_util_dap(gml,1)

	logl 	= J(gml.Nobs[gml.Nlevels,mod],gml.ndim[gml.Nrelevels],0)
	
	//exact events
	gml.survind = 1
	if (merlin_get_nobs(gml)) {
		index 			= merlin_get_surv_index(gml)
		logl[index,] 	= lnnormalden(log(y[index,1]),mu[index,],sigma) :- log(y[index,1])
	}

	//right censored
	gml.survind = 6
	if (merlin_get_nobs(gml)) {
		index 			= merlin_get_surv_index(gml)
		logl[index,] 	= log(1 :- normal((log(y[index,1]):-mu[index,]):/sigma))
	}

	return(logl)
}

`RM' merlin_lognormal_ch(`gml' gml, `RC' t, | `RC' t0)
{
	if (args()==2) 	mu = merlin_util_xzb(gml,t)
	else 			mu = merlin_util_xzb(gml,t,t0)
	sigma	= merlin_util_dap(gml,1)
	
	return(-log(1 :- normal((log(t):-mu):/sigma)))
}

`RM' merlin_lognormal_logh(`gml' gml, `RC' t)
{
	mod		= gml.model
	mu 		= merlin_util_xzb(gml)
	sigma	= merlin_util_dap(gml,1)
	
	logt	= log(t)
	logpdf 	= lnnormalden(logt,mu,sigma) :- logt
	logs 	= log(1 :- normal((logt:-mu):/sigma))
	
	return(logpdf :- logs)
}

end

