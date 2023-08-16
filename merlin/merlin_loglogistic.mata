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

`RM' merlin_logl_loglogistic(`gml' gml, `RM' G, `RM' H)
{
	gml.survind = 0
	model 	= gml.model
	y 	= merlin_util_depvar(gml)
	hasbh	= gml.hasbh[model,]
	Nobs	= merlin_get_nobs(gml,model)
	hast	= gml.istimedep[model,1]
	haslt	= gml.hasltrunc[model]
	logl 	= J(Nobs,1,0)

	mod	= gml.model
	y 	= merlin_util_depvar(gml)
	lambda 	= exp(-merlin_util_xzb(gml))
	gamm	= merlin_util_dap(gml,1)
	onegamm	= 1:/gamm
	
	logl 	= J(merlin_get_nobs(gml,mod),1,0)
	
	//exact events
	gml.survind = 1
	if (merlin_get_nobs(gml)) {
		index 		= merlin_get_surv_index(gml)
		hasbh
		if (hasbh[1] | hasbh[2]) {
			
			logl[index,] = merlin_loglogistic_logh(gml,y)
			
			if (hasbh[1]) {
				logl[index,] = log(exp(logl[index,]) :+ 
					merlin_util_bhazard(gml))
			}
			else {
				logl[index,] = log(exp(logl[index,]) :*
					merlin_util_bhazard(gml))	
			}
		}
		else {
			logl[index,] 	= onegamm:*log(lambda[index,]) :+ 
				(onegamm :- 1):*log(y[index,1]) :- 
				log(gamm) :- 2 :* log(1 :+ 
				(lambda[index,]:*y[index,1]):^onegamm)
		}
	}

	//right censored
	gml.survind = 6
	if (merlin_get_nobs(gml)) {
		index 		= merlin_get_surv_index(gml)
		logl[index,] 	= -log(1 :+ (lambda[index,] :* 
					y[index,1]):^onegamm)
	}

	return(logl)
}

`RM' merlin_logl_loglogistic_ml(`gml' gml)
{
	gml.survind = 0

	mod	= gml.model
	y 	= merlin_util_depvar(gml)
	hasbh	= gml.hasbh[model,]
	lambda 	= exp(-merlin_util_xzb(gml))
	gamm	= merlin_util_dap(gml,1)
	onegamm	= 1:/gamm
	
	logl 	= J(gml.Nobs[gml.Nlevels,mod],gml.ndim[gml.Nrelevels],0)
	
	//exact events
	gml.survind = 1
	if (merlin_get_nobs(gml)) {
		index 		= merlin_get_surv_index(gml)
		if (hasbh[1] | hasbh[2]) {
			
			logl[index,] = merlin_loglogistic_logh(gml,y)
			
			if (hasbh[1]) {
				logl[index,] = log(exp(logl[index,]) :+ 
					merlin_util_bhazard(gml))
			}
			else {
				logl[index,] = log(exp(logl[index,]) :*
					merlin_util_bhazard(gml))	
			}
		} 
		else {
			logl[index,] 	= onegamm:*log(lambda[index,]) :+ 
				(onegamm :- 1):*log(y[index,1]) :- log(gamm) :- 
				2 :* log(1 :+ (lambda[index,] :* 
				y[index,1]):^onegamm)
		}		
	}

	//right censored
	gml.survind = 6
	if (merlin_get_nobs(gml)) {
		index 		= merlin_get_surv_index(gml)
		logl[index,] 	= -log(1 :+ (lambda[index,] :* 
				y[index,1]):^onegamm)
	}

	return(logl)
}

`RM' merlin_loglogistic_ch(`gml' gml, `RC' t, | `RC' t0)
{
	if (args()==2) 	lambda 	= exp(-merlin_util_xzb(gml,t))
	else 			lambda 	= exp(-merlin_util_xzb(gml,t,t0))
	onegamm	= 1:/merlin_util_dap(gml,1)
	
	return(log(1 :+ (lambda:*t):^onegamm))
}

`RM' merlin_loglogistic_logh(`gml' gml, `RC' t)
{
	lambda 	= exp(-merlin_util_xzb(gml))
	gamm	= merlin_util_dap(gml,1)
	onegamm	= 1:/gamm
	
	logpdf 	= onegamm:*log(lambda) :+ (onegamm :- 1) :* 
			log(t) :- log(gamm) :- 2 :* 
			log(1 :+ (lambda:*t):^onegamm)
	logs 	= -log(1 :+ (lambda:*t):^onegamm)
	
	return(logpdf :- logs)
}

end

