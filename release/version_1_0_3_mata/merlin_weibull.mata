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

`RM' merlin_logl_weibull(`gml' gml)
{
	y 		= merlin_util_depvar(gml)
	
	xzb	 	= merlin_util_xzb(gml,y[,1])
	gam 	= asarray(gml.distancb,(gml.model,1))

	logh 	= merlin_weibull_logh(xzb, gam, y[,1])
	
	if (gml.nobhaz[gml.model]) {
		if (gml.NI[gml.model]) 	result = y[,2] :* logh :- merlin_weibull_ch_ni(gml,gam,y[,1])
		else 					result = y[,2] :* logh :- merlin_weibull_ch(xzb,gam,y[,1]) 
	}
	else {
		if (gml.NI[gml.model]) {
			result = y[,2] :* log(exp(logh) :+ asarray(gml.bhazards,gml.model)) :- merlin_weibull_ch_ni(gml,gam,y[,1])
		}
		else {
			result = y[,2] :* log(exp(logh) :+ asarray(gml.bhazards,gml.model)) :- merlin_weibull_ch(xzb,gam,y[,1]) 
		}
	}

	return(result)
}

`RM' merlin_weibull_expval(`gml' gml, | `RC' t)
{
	gam 	= asarray(gml.distancb,(gml.model,1))
	if (args()==1) 	return(gamma(1+1/gam) :* exp(-merlin_util_xzb(gml):/gam))
	else			return(gamma(1+1/gam) :* exp(-merlin_util_xzb(gml,t):/gam))
}

`RM' merlin_weibull_logh(`RM' xzb, `RS' gam, `RC' t)
{
	return(xzb :+ log(gam) :+ (gam - 1) :* log(t))
}

`RM' merlin_weibull_logh_off(`RM' xzb, `RS' gam, `RC' t)
{
	return(xzb :+ log(gam) :+ gam :* log(t))
}

`RM' merlin_weibull_ch(`RM' xzb, `RS' gam, `RC' t)
{
	return(exp(xzb) :* t :^ gam)
}

`RM' merlin_weibull_ch_ni(`gml' gml, `RS' gam, | `RC' t)
{
	cumhaz 		= J(gml.Nobs[gml.Nlevels,gml.model],1,0)
	hq			= asarray(gml.haznodes,gml.model)
	loghazq 	= log(gam) :+ (gam-1) :* log(hq) :+ log(asarray(gml.hazweights,gml.model))
	
	for (q=1;q<=gml.hazNnodes[gml.model];q++) {
		cumhaz = cumhaz :+ exp(merlin_util_xzb(gml,hq[,q]) :+ loghazq[,q])
	}
	
	return(cumhaz)
}

//for left truncation
`RM' merlin_weibull_s(`gml' gml, | `RC' t)
{
	
	gam = asarray(gml.distancb,(gml.model,1))
	if (gml.NI[gml.model])	return(exp(-merlin_weibull_ch_ni(gml,gam,t))) 
	else 					return(exp(-merlin_weibull_ch(gml,gam,t))) 
	
}

end
