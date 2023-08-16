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

`RM' merlin_logl_ggamma(`gml' gml, `RM' G, `RM' H)
{
	gml.survind = 0

	mod	= gml.model
	y 	= merlin_util_depvar(gml)
	mu 	= merlin_util_xzb(gml)
	sigma	= merlin_util_dap(gml,1)
	kapp	= merlin_util_dap(gml,2)
	gamm	= abs(kapp)^(-2)
	
	z	= sign(kapp) :* (log(y[,1]) :- mu):/sigma
	u	= gamm :* exp(abs(kapp) :* z)
	
	logl 	= J(merlin_get_nobs(gml,mod),1,0)
	
	//exact events
	gml.survind = 1
	if (merlin_get_nobs(gml)) {
		index = merlin_get_surv_index(gml)
		if (kapp==0) {
			logl[index,] = -log(sigma :* y[index,1] :* sqrt(2*pi())) :- (z[index,]:^2):/2
		}
		else {
			logl[index,] = gamm:*log(gamm) :- log(sigma :* y[index,1] :* sqrt(gamm) :* gamma(gamm)) :+ z[index,]:*sqrt(gamm) :- u[index,]
		}
	}

	//right censored
	gml.survind = 6
	if (merlin_get_nobs(gml)) {
		index = merlin_get_surv_index(gml)
		if (kapp>0) {
			logl[index,] = log(1 :- gammap(gamm,u[index,]))
		}
		else if (kapp<0) {
			logl[index,] = log(gammap(gamm,u[index,]))
		}
		else {
			logl[index,] = log(1 :- normal(z[index,]))
		}
	}

	return(logl)
}


`RM' merlin_logl_ggamma_ml(`gml' gml)
{
	gml.survind = 0

	mod		= gml.model
	y 		= merlin_util_depvar(gml)
	mu 		= merlin_util_xzb(gml)
	sigma	= merlin_util_dap(gml,1)
	kapp	= merlin_util_dap(gml,2)
	gamm	= abs(kapp)^(-2)
	
	z		= sign(kapp) :* (log(y[,1]) :- mu):/sigma
	u		= gamm :* exp(abs(kapp) :* z)
	
	logl 	= J(gml.Nobs[gml.Nlevels,mod],gml.ndim[gml.Nrelevels],0)
	
	//exact events
	gml.survind = 1
	if (merlin_get_nobs(gml)) {
		index = merlin_get_surv_index(gml)
		if (kapp==0) {
			logl[index,] = -log(sigma :* y[index,1] :* sqrt(2*pi())) :- (z[index,]:^2):/2
		}
		else {
			logl[index,] = gamm:*log(gamm) :- log(sigma :* y[index,1] :* sqrt(gamm) :* gamma(gamm)) :+ z[index,]:*sqrt(gamm) :- u[index,]
		}
	}

	//right censored
	gml.survind = 6
	if (merlin_get_nobs(gml)) {
		index = merlin_get_surv_index(gml)
		if (kapp>0) {
			logl[index,] = log(1 :- gammap(gamm,u[index,]))
		}
		else if (kapp<0) {
			logl[index,] = log(gammap(gamm,u[index,]))
		}
		else {
			logl[index,] = log(1 :- normal(z[index,]))
		}
	}

	return(logl)
}

`RM' merlin_ggamma_ch(`gml' gml, `RC' t, | `RC' t0)
{
	mod		= gml.model
	if (args()==2) 	mu = merlin_util_xzb(gml,t)
	else			mu = merlin_util_xzb(gml,t,t0)
	sigma	= merlin_util_dap(gml,1)
	kapp	= merlin_util_dap(gml,2)
	gamm	= abs(kapp)^(-2)
	z		= sign(kapp) :* (log(t) :- mu):/sigma
	u		= gamm :* exp(abs(kapp) :* z)
	
	if (kapp>0) 		return(-log(1 :- gammap(gamm,u)))
	else if (kapp<0) 	return(-log(gammap(gamm,u)))
	else 				return(-log(1 :- normal(z)))
}

`RM' merlin_ggamma_logh(`gml' gml, `RC' t)
{
	mod		= gml.model
	y 		= merlin_util_depvar(gml)
	mu 		= merlin_util_xzb(gml)
	sigma	= merlin_util_dap(gml,1)
	kapp	= merlin_util_dap(gml,2)
	gamm	= abs(kapp)^(-2)	
	z		= sign(kapp) :* (log(t) :- mu):/sigma
	u		= gamm :* exp(abs(kapp) :* z)
	
	if (kapp==0) logpdf = -log(sigma :* t :* sqrt(2*pi())) :- (z:^2):/2
	else		 logpdf = gamm:*log(gamm) :- log(sigma :* t :* sqrt(gamm) :* gamma(gamm)) :+ z:*sqrt(gamm) :- u

	if (kapp>0) 		logs = log(1 :- gammap(gamm,u))
	else if (kapp<0) 	logs = log(gammap(gamm,u))
	else 				logs = log(1 :- normal(z))
	
	return(logpdf :- logs)
}

end

