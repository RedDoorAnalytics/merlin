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

`RM' merlin_logl_exp(`gml' gml)
{
	y 		= merlin_util_depvar(gml)
	xzb		= merlin_util_xzb(gml,y[,1])
	logh	= xzb
	
	if (gml.nobhaz[gml.model]) {
		if (gml.NI[gml.model]) 	{
			result = y[,2] :* logh :- merlin_exp_ch_ni(gml,y[,1])
		}
		else {
			result = y[,2] :* logh :- merlin_exp_ch(xzb,y[,1]) 
		}
	}
	else {
		if (gml.NI[gml.model]) {
			result = y[,2] :* log(exp(logh) :+ asarray(gml.bhazards,gml.model)) :- merlin_exp_ch_ni(gml,y[,1])
		}
		else {
			result = y[,2] :* log(exp(logh) :+ asarray(gml.bhazards,gml.model)) :- merlin_exp_ch(xzb,y[,1]) 
		}
	}
	return(result)
}

`RM' merlin_exp_expval(`gml' gml, | `RC' t)
{
	if (args()==1) 	return(1:/exp(merlin_util_xzb(gml)))
	else	 		return(1:/exp(merlin_util_xzb(gml,t)))
}

`RM' merlin_exp_h(`RM' xzb, `RS' gam, `RC' t)
{
	return(exp(xzb))
}

`RM' merlin_exp_logh(`RM' xzb, `RC' t)
{
	return(xzb)
}

`RM' merlin_exp_logh_off(`RM' xzb,`RC' t)
{
	return(xzb :+ log(t))
}

`RM' merlin_exp_ch(`RM' xzb, `RC' t)
{
	return(exp(xzb) :* t)
}

`RM' merlin_exp_ch_ni(`gml' gml, | `RC' t)
{
	cumhaz 		= J(gml.Nobs[gml.Nlevels,gml.model],1,0)
	hq			= asarray(gml.haznodes,gml.model)
	loghazq 	= log(asarray(gml.hazweights,gml.model))
	
	for (q=1;q<=gml.hazNnodes[gml.model];q++) {
		cumhaz = cumhaz :+ exp(merlin_util_xzb(gml,hq[,q]) :+ loghazq[,q])
	}
	
	return(cumhaz)
}

`RM' merlin_exp_s(`RM' xzb, `RC' t)
{
	return(exp(-exp(xzb) :* t ))
}

`RM' merlin_exp_s_ni(`gml' gml, `RC' t)
{
	return(exp(-merlin_exp_ch_ni(gml,t)))
}

end
