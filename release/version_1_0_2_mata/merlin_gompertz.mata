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

`RM' merlin_logl_gompertz(`gml' gml)
{
	y 		= merlin_util_depvar(gml)
	
	xzb	 	= merlin_util_xzb(gml,y[,1])
	gam 	= asarray(gml.distancb,(gml.model,1))
	
	logh 	= merlin_gompertz_logh(xzb, gam, y[,1])
	
	if (gml.nobhaz[gml.model]) {
		if (gml.NI[gml.model]) 	{
			result = y[,2] :* logh :- merlin_gompertz_ch_ni(gml,gam,y[,1])
		}
		else {
			result = y[,2] :* logh :- merlin_gompertz_ch(xzb,gam,y[,1]) 
			if (gml.hasltrunc[gml.model]) result = result :+ merlin_gompertz_ch(merlin_util_xzb(gml,y[,3]),gam,y[,3])
		}
	}
	else {
		if (gml.NI[gml.model]) {
			result = y[,2] :* log(exp(logh) :+ asarray(gml.bhazards,gml.model)) :- merlin_gompertz_ch_ni(gml,gam,y[,1])
		}
		else {
			result = y[,2] :* log(exp(logh) :+ asarray(gml.bhazards,gml.model)) :- merlin_gompertz_ch(xzb,gam,y[,1]) 
			if (gml.hasltrunc[gml.model]) result = result :+ merlin_gompertz_ch(merlin_util_xzb(gml,y[,3]),gam,y[,1])
		}
	}
	return(result)
}

`RM' merlin_gompertz_h(`RM' xzb, `RS' gam, `RC' t)
{
	return(exp(xzb :+ gam :* t))
}

`RM' merlin_gompertz_logh(`RM' xzb, `RS' gam, `RC' t)
{
	return(xzb :+ gam :* t)
}

`RM' merlin_gompertz_logh_off(`RM' xzb, `RS' gam, `RC' t)
{
	return(xzb :+ gam :* t :+ log(t))
}

`RM' merlin_gompertz_ch(`RM' xzb, `RS' gam, `RC' t)
{
	return(exp(xzb) :* (1/gam) :*(exp(gam:*t):-1))
}

`RM' merlin_gompertz_ch_ni(`gml' gml, `RS' gam, | `RC' t)
{
	cumhaz 		= J(gml.Nobs[gml.Nlevels,gml.model],1,0)
	hq			= asarray(gml.haznodes,gml.model)
	loghazq 	= gam :* hq :+ log(asarray(gml.hazweights,gml.model))
	
	for (q=1;q<=gml.hazNnodes[gml.model];q++) {
		cumhaz = cumhaz :+ exp(merlin_util_xzb(gml,hq[,q]) :+ loghazq[,q])
	}
	
	return(cumhaz)
}

`RM' merlin_gompertz_s(`RM' xzb, `RS' gam, `RC' t)
{
	return(exp(-merlin_gompertz_ch(xzb,gam,t)))
}

`RM' merlin_gompertz_s_ni(`gml' gml, `RS' gam, `RC' t)
{
	return(exp(-merlin_gompertz_ch_ni(gml,gam,t)))
}

end
